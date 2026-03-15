#! /usr/bin/env -S python3 -i
"""
EPFRAME Visualization Tools
Plots frame geometry, deformed shape, and force diagrams

Reads all data from the output file (which contains input echo)
Uses 0-indexed numpy arrays internally, 1-indexed for display
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, FancyBboxPatch, FancyArrowPatch
import sys
import re
import os

def read_output_file(filename):
    """
    Read output file and extract all data including input echo and results.
    Handles both old numeric DOF format (0/1) and new string format (BI/POS/NEG/FREE).
    Returns geometry, element data, RTYPE, and analysis results including lift-off history.
    """
    _RTYPE_MAP = {'FREE': 0, 'BI': 1, 'POS': 2, 'NEG': 3}

    with open(filename, 'r') as f:
        lines = f.readlines()

    FN = 0        # kept for backward compat; new output files have no frame number
    title = ''    # free-text title from first header line
    NCT = 0
    NE = 0
    E_mod = 0.0
    Fy    = 0.0
    DI    = None  # degree of static indeterminacy (None if old output format)
    DLMT  = None  # displacement limit (None if old output format)
    CORD = None
    ECON = None
    PM = None
    SMA = None
    DOF = None
    RTYPE = None

    # --- Pass 1: parse header / input echo ---
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # New format: title is the first non-blank % line after the opening %
        # Old format: contains 'FRAME NO'
        if 'FRAME NO' in line:
            parts = line.split()
            for j, p in enumerate(parts):
                if p == 'NO':
                    FN = int(parts[j+1])
                    break
        # Title line: a % line that is NOT a separator or known keyword
        if (line.startswith('%') and title == '' and
                'ELASTIC' not in line and 'FRAME' not in line and
                'DATA' not in line and '---' not in line and
                len(line.replace('%', '').strip()) > 0):
            title = line.replace('%', '').strip()

        if 'NUMBER OF NODES'      in line: NCT   = int(line.split()[-1])
        if 'NUMBER OF ELEMENTS'   in line: NE    = int(line.split()[-1])
        if 'MOD OF ELASTICITY'    in line: E_mod = float(line.split()[-1])
        if 'YIELD STRESS'         in line: Fy    = float(line.split()[-1])
        # New header lines: "STATIC INDETERMINACY   3  (mechanism..." → value is 3rd token
        if 'STATIC INDETERMINACY' in line:
            try: DI = int(line.replace('%','').split()[2])
            except (ValueError, IndexError): pass
        if 'DISPLACEMENT LIMIT'   in line:
            try: DLMT = float(line.replace('%','').split()[2])
            except (ValueError, IndexError): pass

        if 'DATA FOR NODES' in line:
            i += 2  # skip section header and column-header line
            CORD  = np.zeros((NCT, 2))
            DOF   = np.zeros((NCT, 3), dtype=int)
            RTYPE = np.zeros((NCT, 3), dtype=int)

            for nd in range(NCT):
                i += 1
                data_line = lines[i].replace('%', '').strip()
                parts = data_line.split()
                node_num = int(parts[0]) - 1  # 0-indexed

                CORD[node_num, 0] = float(parts[1])
                CORD[node_num, 1] = float(parts[2])

                # Support new string format (BI/POS/NEG/FREE) and old integer format (0/1)
                for coord in range(3):
                    token = parts[3 + coord]
                    if token in _RTYPE_MAP:
                        # New format
                        rt = _RTYPE_MAP[token]
                    else:
                        # Old format: 0 = free, 1 = bidirectional reaction
                        rt = int(token)   # 1 → BI, 0 → FREE
                    RTYPE[node_num, coord] = rt
                    # Bidirectional supports are permanently constrained (excluded from DOF vec)
                    # Unidirectional supports remain in DOF vec (active-set handles them)
                    DOF[node_num, coord] = 0 if rt == 1 else 1

        if 'YIELD STRESS'      in line: Fy    = float(line.split()[-1])

        if 'DATA FOR ELEMENTS' in line:
            i += 2
            ECON = np.zeros((NE, 2), dtype=int)
            PM   = np.zeros(NE)
            SMA  = np.zeros(NE)   # moment of inertia IXX per element

            for el in range(NE):
                i += 1
                data_line = lines[i].replace('%', '').strip()
                parts = data_line.split()
                elem_num = int(parts[0]) - 1
                ECON[elem_num, 0] = int(parts[1]) - 1
                ECON[elem_num, 1] = int(parts[2]) - 1
                SMA[elem_num]     = float(parts[3])   # IXX
                # New format: parts[5]=Z, parts[6]=MP  (7-column echo)
                # Old format: parts[5]=MP              (6-column echo)
                if len(parts) >= 7:
                    PM[elem_num] = float(parts[6])    # MP = Z * Fy
                else:
                    PM[elem_num] = float(parts[5])    # backward compat

        if 'PLASTIC HINGE' in line or (not line.startswith('%') and line and 'NCYCL' not in line):
            break

        i += 1

    # --- Pass 2: parse analysis results ---
    hinges        = []
    disp_history  = []
    moment_history = []
    force_history  = []
    load_factors  = []
    liftoff_history = []   # list of sets: {(node_0idx, coord_name), ...} per stage

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if 'PLASTIC HINGE' in line and 'FORMED' in line:
            parts = line.replace('%', '').split()

            try:
                hinge_num   = int(parts[parts.index('HINGE')   + 1])
                elem_num    = int(parts[parts.index('ELEMENT')  + 1]) - 1
                node_num    = int(parts[parts.index('NODE')     + 1]) - 1
                load_factor = float(parts[parts.index('IS')     + 1])
            except (ValueError, IndexError):
                print(f"Warning: Could not parse hinge line: {line}")
                i += 1
                continue

            hinges.append((hinge_num, elem_num, node_num))
            load_factors.append(load_factor)

            # --- active support status (lift-off) ---
            liftoff_nodes = set()
            j = i + 1
            while j < len(lines):
                sl = lines[j].strip()
                if 'ACTIVE SUPPORT STATUS' in sl:
                    j += 1
                    while j < len(lines):
                        sl2 = lines[j].strip().replace('%', '').strip()
                        # Lines look like:  NODE   4: Y=POS:LIFT-OFF
                        if not sl2 or 'CUMULATIVE' in sl2 or 'REACTIONS' in sl2:
                            break
                        if 'NODE' in sl2 and 'LIFT-OFF' in sl2:
                            m = re.search(r'NODE\s+(\d+):(.*)', sl2)
                            if m:
                                nd_1 = int(m.group(1)) - 1  # 0-indexed
                                # parse individual DOF tokens  e.g. Y=POS:LIFT-OFF
                                for token in m.group(2).split(','):
                                    if 'LIFT-OFF' in token:
                                        coord = token.strip().split('=')[0]
                                        liftoff_nodes.add((nd_1, coord))
                        j += 1
                    break
                if 'CUMULATIVE' in sl:
                    break
                j += 1
            liftoff_history.append(liftoff_nodes)

            # --- cumulative deformations ---
            while i < len(lines) and 'X-DISP' not in lines[i]:
                i += 1
            i += 1

            disp = np.zeros((NCT, 3))
            for j in range(NCT):
                if i >= len(lines): break
                data_line = lines[i].strip().replace('%', '')
                if not data_line or 'CUMULATIVE' in data_line or 'ELEMENT' in data_line:
                    break
                parts = data_line.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    nd = int(parts[0]) - 1
                    disp[nd, 0] = float(parts[1])
                    disp[nd, 1] = float(parts[2])
                    disp[nd, 2] = float(parts[3])
                i += 1
            disp_history.append(disp)

            # --- cumulative moments ---
            while i < len(lines) and 'END MOMENTS' not in lines[i]:
                i += 1
            i += 1

            moments = np.zeros((NE, 2))
            for j in range(NE):
                if i >= len(lines): break
                data_line = lines[i].strip().replace('%', '')
                if not data_line or 'CUMULATIVE' in data_line or 'TENSION' in data_line:
                    break
                parts = data_line.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    el = int(parts[0]) - 1
                    moments[el, 0] = float(parts[1])
                    moments[el, 1] = float(parts[2])
                i += 1
            moment_history.append(moments)

            # --- cumulative tension forces ---
            while i < len(lines) and 'CUMULATIVE TENSION' not in lines[i]:
                i += 1
            i += 2   # skip "CUMULATIVE TENSION FORCES" and "ELEMENT  TENSION" header

            forces = np.zeros(NE)
            for j in range(NE):
                if i >= len(lines): break
                data_line = lines[i].strip().replace('%', '')
                if not data_line or 'ANALYSIS' in data_line or 'PLASTIC' in data_line \
                        or 'REACTIONS' in data_line:
                    break
                parts = data_line.split()
                if len(parts) >= 2 and parts[0].isdigit():
                    el = int(parts[0]) - 1
                    forces[el] = float(parts[1])
                i += 1
            force_history.append(forces)

        i += 1

    return (FN, NCT, NE, E_mod, CORD, DOF, RTYPE, ECON, PM, SMA,
            hinges, disp_history, moment_history, force_history,
            load_factors, liftoff_history, title, DI, DLMT)

def compute_element_curves(CORD, ECON, SMA, E_mod, disp, moments, NE,
                           frame_size, scale=1.0):
    """
    Compute the curved deformed shape for each element using cubic Hermite
    interpolation driven by the slope-deflection flexibility relation.

    For point loads at nodes only the moment varies linearly along each element
    and the elastic curve is a cubic.  Element end rotations relative to the
    deformed chord follow from:

        φ₁ = (L / 6EI)(2M₁ − M₂)   [near-end, node n1]
        φ₂ = (L / 6EI)(2M₂ − M₁)   [far-end,  node n2]

    At a plastic hinge the computed φ differs from the structural node rotation
    in CD, so the rotation discontinuity appears automatically in the plotted
    curve — no special hinge logic is required.

    INPUTS
        CORD        (NCT,2)  original node coordinates
        ECON        (NE,2)   element connectivity (0-indexed)
        SMA         (NE,)    second moment of area IXX per element
        E_mod                modulus of elasticity
        disp        (NCT,3)  cumulative nodal displacements [ΔX, ΔY, Δθ] — UNSCALED
        moments     (NE,2)   cumulative end moments [M_near, M_far]
        NE                   number of elements
        frame_size           max(x_span, y_span) of the undeformed frame — sets dx
        scale                displacement amplification factor for visualisation

    OUTPUTS
        curves  list of NE arrays, each shape (n_pts, 2) giving global (x, y)
                of the curved deformed element sampled at dx ≈ 0.01 × frame_size
    """
    dx_step = 0.01 * frame_size   # target arc-length step in original coordinates

    curves = []

    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]

        # Scaled deformed end positions (chord endpoints for the visualisation)
        X1 = CORD[n1, 0] + scale * disp[n1, 0]
        Y1 = CORD[n1, 1] + scale * disp[n1, 1]
        X2 = CORD[n2, 0] + scale * disp[n2, 0]
        Y2 = CORD[n2, 1] + scale * disp[n2, 1]

        # Original element length (used for EI and rotation calculation)
        dx0 = CORD[n2, 0] - CORD[n1, 0]
        dy0 = CORD[n2, 1] - CORD[n1, 1]
        L0  = np.sqrt(dx0**2 + dy0**2)

        # Deformed chord vector
        dX = X2 - X1
        dY = Y2 - Y1
        L_def = np.sqrt(dX**2 + dY**2)

        # Number of segments: ceil(L0 / dx_step), minimum 2
        n_segs = max(2, int(np.ceil(L0 / dx_step)))
        t_vals = np.linspace(0.0, 1.0, n_segs + 1)

        if L_def < 1.0e-12 or L0 < 1.0e-12 or SMA[el] < 1.0e-30:
            curves.append(np.column_stack([
                X1 + t_vals * dX,
                Y1 + t_vals * dY]))
            continue

        # Unit perpendicular to deformed chord (90° CCW)
        nx = -dY / L_def
        ny =  dX / L_def

        # Hermite basis functions for relative transverse displacement
        H2 = t_vals * (1.0 - t_vals)**2   # near-end rotation contribution
        H4 = t_vals**2 * (t_vals - 1.0)   # far-end  rotation contribution

        # Slope-deflection flexibility → relative end rotations (physical)
        EI   = E_mod * SMA[el]
        M1, M2 = moments[el, 0], moments[el, 1]
        fac  = L0 / (6.0 * EI)
        phi1 = fac * (2.0*M1 - M2)   # near-end (node n1)
        phi2 = fac * (2.0*M2 - M1)   # far-end  (node n2)

        # Transverse displacement relative to deformed chord — scaled for visualisation
        v_rel = scale * L0 * (H2 * phi1 + H4 * phi2)

        # Global coordinates along the curved element
        x_curve = X1 + t_vals * dX + v_rel * nx
        y_curve = Y1 + t_vals * dY + v_rel * ny

        curves.append(np.column_stack([x_curve, y_curve]))

    return curves


def _rxn_arrow(ax, x0, y0, x1, y1, both_ends=False, zorder=6):
    """
    Draw a reaction line using FancyArrowPatch.
    shrinkA=shrinkB=0 ensures the tip reaches the exact endpoint in data coords.
    mutation_scale sets the arrowhead size in display points.
    both_ends=True adds a head at both tips (bidirectional reaction).
    """
    kw = dict(arrowstyle='->', color='red', linewidth=2.0,
              mutation_scale=20, shrinkA=0, shrinkB=0, zorder=zorder)
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), **kw))
    if both_ends:
        ax.add_patch(FancyArrowPatch((x1, y1), (x0, y0), **kw))


def draw_support(ax, x, y, RTYPE, d_node):
    """
    Draw reaction support icon centred on the node.
    Used on both the geometry plot and the deformed-shape plot.

    Every icon has:
      • Base shape: solid red circle (no Z reaction) or solid red square (Z reaction)
      • Horizontal line of total length 4×d_node centred at the node
      • Vertical   line of total length 4×d_node centred at the node

    Arrowhead convention on lines:
      Reaction present (rtype > 0): FancyArrowPatch with head(s) at tip(s)
        rtype 1 (bidirectional) → heads at both ends
        rtype 2 (+)             → head at positive end
        rtype 3 (−)             → head at negative end
      No reaction in that direction (rtype == 0): plain red line, no arrowhead

    Sizes:
      d_node : node diameter in data coords (caller sets margin × 0.08)
      d = 3 × d_node : icon diameter / side
      r = d / 2
      hl = 2 × d_node : half-length of lines  (total = 4 × d_node)
    """
    rx, ry, rz = int(RTYPE[0]), int(RTYPE[1]), int(RTYPE[2])
    has_x = rx > 0;  has_y = ry > 0;  has_z = rz > 0

    d      = 2.0 * d_node
    r      = d / 2.0
    hl     = 3.0 * d_node   # half-length: tips at 3×d_node, icon edge at 1×d_node
    RED    = 'red'
    Z_ICON = 5
    Z_LINE = 6              # lines/arrows drawn above the icon patch

    # --- base icon at Z_ICON ---
    if has_z:
        sq = Polygon([(x-r, y-r), (x+r, y-r), (x+r, y+r), (x-r, y+r)],
                     closed=True, facecolor=RED, edgecolor=RED,
                     linewidth=1.5, zorder=Z_ICON)
        ax.add_patch(sq)
    else:
        circ = plt.Circle((x, y), r, facecolor=RED, edgecolor=RED,
                           linewidth=1.5, zorder=Z_ICON)
        ax.add_patch(circ)

    # --- horizontal line at Z_LINE (above icon) ---
    x_L, x_R = x - hl, x + hl
    if rx == 1:
        _rxn_arrow(ax, x_L, y, x_R, y, both_ends=True, zorder=Z_LINE)
    elif rx == 2:
        _rxn_arrow(ax, x_L, y, x_R, y, zorder=Z_LINE)
    elif rx == 3:
        _rxn_arrow(ax, x_R, y, x_L, y, zorder=Z_LINE)
    else:
        ax.plot([x_L, x_R], [y, y], color=RED, linewidth=2.0, zorder=Z_LINE)

    # --- vertical line at Z_LINE (above icon) ---
    y_B, y_T = y - hl, y + hl
    if ry == 1:
        _rxn_arrow(ax, x, y_B, x, y_T, both_ends=True, zorder=Z_LINE)
    elif ry == 2:
        _rxn_arrow(ax, x, y_B, x, y_T, zorder=Z_LINE)
    elif ry == 3:
        _rxn_arrow(ax, x, y_T, x, y_B, zorder=Z_LINE)
    else:
        ax.plot([x, x], [y_B, y_T], color=RED, linewidth=2.0, zorder=Z_LINE)


# Identical icon used on both geometry and deformed-shape plots
draw_reaction_arrows = draw_support



def plot_frame_geometry(CORD, DOF, RTYPE, ECON, NCT, NE, title="Frame Geometry"):
    """Plot undeformed frame geometry with supports"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Find plot limits
    x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
    y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    margin = max(x_range, y_range) * 0.15
    
    # Draw elements
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        x = [CORD[n1, 0], CORD[n2, 0]]
        y = [CORD[n1, 1], CORD[n2, 1]]
        ax.plot(x, y, 'b-', linewidth=3, zorder=1)
        
        # Label element at midpoint
        mid_x = (CORD[n1, 0] + CORD[n2, 0]) / 2
        mid_y = (CORD[n1, 1] + CORD[n2, 1]) / 2
        ax.text(mid_x, mid_y, f'{el+1}', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7),
                ha='center', va='center', zorder=2)
    
    # Node reference diameter in data coordinates
    d_node = margin * 0.08

    # Draw reaction icons (before nodes so nodes overlap them)
    for nd in range(NCT):
        xn, yn = CORD[nd, 0], CORD[nd, 1]
        if np.any(RTYPE[nd] > 0):
            draw_support(ax, xn, yn, RTYPE[nd], d_node)

    # Draw nodes last so they sit on top of everything (zorder=7)
    for nd in range(NCT):
        xn, yn = CORD[nd, 0], CORD[nd, 1]
        ax.plot(xn, yn, 'ko', markersize=8, zorder=7)
        ax.text(xn + margin * 0.07, yn + margin * 0.07, f'{nd+1}',
                fontsize=11, fontweight='bold',
                ha='left', va='bottom', color='black', zorder=8)
    
    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (in)', fontsize=12)
    ax.set_ylabel('Y (in)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    return fig, ax

def plot_deformed_shape(CORD, DOF, RTYPE, ECON, NCT, NE, disp, moments, load_factor,
                        SMA, E_mod, scale=None, hinge_info=None, liftoff_nodes=None):
    """Plot deformed shape with optional hinge markers"""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Find plot limits
    x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
    y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    margin  = max(x_range, y_range) * 0.15
    frame_size = max(x_range, y_range)
    d_node  = margin * 0.08

    # Auto-scale based on max nodal translation
    if scale is None:
        max_disp = np.max(np.abs(disp[:, :2]))
        scale = frame_size * 0.08 / max_disp if max_disp > 0 else 1.0

    # Draw original shape (dashed)
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        ax.plot([CORD[n1, 0], CORD[n2, 0]], [CORD[n1, 1], CORD[n2, 1]],
                'b--', linewidth=1.5, alpha=0.4, label='Original' if el == 0 else '')

    # Draw curved deformed shape
    curves = compute_element_curves(CORD, ECON, SMA, E_mod,
                                    disp, moments, NE, frame_size, scale=scale)
    for el, curve in enumerate(curves):
        ax.plot(curve[:, 0], curve[:, 1], 'r-', linewidth=2.5,
                label='Deformed' if el == 0 else '')

    # Reaction icons at original support positions (same as geometry plot)
    for nd in range(NCT):
        if np.any(RTYPE[nd] > 0):
            draw_support(ax, CORD[nd, 0], CORD[nd, 1], RTYPE[nd], d_node)

    # Draw plastic hinge markers on the deformed curves
    if hinge_info:
        for h_num, el_num, nd_num in hinge_info:
            curve = curves[el_num]
            if nd_num == ECON[el_num, 0]:
                x_h, y_h = curve[0, 0], curve[0, 1]
            else:
                x_h, y_h = curve[-1, 0], curve[-1, 1]
            circle = Circle((x_h, y_h), margin * 0.06,
                            facecolor='white', edgecolor='red', linewidth=2.5, zorder=10)
            ax.add_patch(circle)
            ax.text(x_h, y_h, f'{h_num}', fontsize=9, fontweight='bold',
                    color='red', ha='center', va='center', zorder=11)

    # Draw deformed node positions — label upper-right to avoid overlap with reaction lines
    for nd in range(NCT):
        x_def = CORD[nd, 0] + scale * disp[nd, 0]
        y_def = CORD[nd, 1] + scale * disp[nd, 1]
        ax.plot(x_def, y_def, 'ro', markersize=6, zorder=7)
        ax.text(x_def + margin * 0.07, y_def + margin * 0.07, f'{nd+1}',
                fontsize=10, fontweight='bold',
                ha='left', va='bottom', color='black', zorder=8)
    
    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(r'$X$ (in)', fontsize=12)
    ax.set_ylabel(r'$Y$ (in)', fontsize=12)
    ax.set_title(rf'Deformed Shape ($\lambda$ = {load_factor:.3f}, scale = {scale:.1f}$\times$)',
                fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    
    # Displacement summary
    max_x = np.max(np.abs(disp[:, 0]))
    max_y = np.max(np.abs(disp[:, 1]))
    max_r = np.max(np.abs(disp[:, 2]))
    info_text = (f'Max displacements:\n'
                 rf'$\Delta X$: {max_x:.4f} in'  '\n'
                 rf'$\Delta Y$: {max_y:.4f} in'  '\n'
                 rf'$\theta$: {max_r:.5f} rad')
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    return fig, ax

def plot_moment_diagram(CORD, ECON, NCT, NE, moments, PM, load_factor, hinges=None):
    """Plot bending moment diagram"""
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Find plot limits
    x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
    y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    margin = max(x_range, y_range) * 0.15
    
    # Find max moment for scaling
    max_moment = np.max(np.abs(moments))
    if max_moment > 0:
        moment_scale = margin * 0.8 / max_moment
    else:
        moment_scale = 1.0
    
    # Draw elements and moment diagrams
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        x1, y1 = CORD[n1, 0], CORD[n1, 1]
        x2, y2 = CORD[n2, 0], CORD[n2, 1]

        # Draw element
        ax.plot([x1, x2], [y1, y2], 'b-', linewidth=2.5, zorder=1)

        Mi = moments[el, 0]   # element-coord near-end moment (CCW positive)
        Mj = moments[el, 1]   # element-coord far-end moment

        # Internal moment: M(x) = -Mi*(1-s) + Mj*s
        # Positive internal moment → compression on top → offset in 90°-CCW normal direction
        dx = x2 - x1;  dy = y2 - y1
        L  = np.sqrt(dx**2 + dy**2)
        nx = -dy / L;  ny = dx / L   # 90° CCW = left of element direction

        # Moment diagram with correct sign at near end
        n_points = 40
        x_diag, y_diag = [], []
        for j in range(n_points + 1):
            s    = j / n_points
            M_pt = -Mi * (1.0 - s) + Mj * s   # internal moment
            x_diag.append(x1 + s*dx + M_pt * moment_scale * nx)
            y_diag.append(y1 + s*dy + M_pt * moment_scale * ny)

        x_full = [x1] + x_diag + [x2]
        y_full = [y1] + y_diag + [y2]
        ax.fill(x_full, y_full, color='lightgray', alpha=0.6, edgecolor='blue', linewidth=1.5)

        # End-value labels: internal moment at near end = -Mi, at far end = Mj
        M_near = -Mi
        M_far  =  Mj
        lbl_off = 1.2
        if abs(M_near) > 1e-3:
            ax.text(x1 + M_near*moment_scale*nx*lbl_off,
                    y1 + M_near*moment_scale*ny*lbl_off,
                    f'{M_near:.0f}', fontsize=9, ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        if abs(M_far) > 1e-3:
            ax.text(x2 + M_far*moment_scale*nx*lbl_off,
                    y2 + M_far*moment_scale*ny*lbl_off,
                    f'{M_far:.0f}', fontsize=9, ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        # Mark plastic hinges
        if hinges:
            for h_num, el_num, nd_num in hinges:
                if el_num == el:
                    x_h, y_h = (x1, y1) if nd_num == n1 else (x2, y2)
                    circle = Circle((x_h, y_h), margin * 0.06,
                                    facecolor='white', edgecolor='red',
                                    linewidth=2.5, zorder=10)
                    ax.add_patch(circle)
                    ax.text(x_h, y_h - margin * 0.12, f'H{h_num}',
                            fontsize=9, fontweight='bold', color='red',
                            ha='center', va='top')

    # Draw nodes (no labels on force diagrams)
    for nd in range(NCT):
        ax.plot(CORD[nd, 0], CORD[nd, 1], 'ko', markersize=8, zorder=3)

    ax.set_xlim(x_min - margin * 1.5, x_max + margin * 1.5)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin * 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(r'$X$ (in)', fontsize=12)
    ax.set_ylabel(r'$Y$ (in)', fontsize=12)
    ax.set_title(rf'Bending Moment Diagram ($\lambda$ = {load_factor:.3f})',
                 fontsize=14, fontweight='bold')

    # Legend — no plastic moment value
    ax.text(0.02, 0.98,
            f'Peak |M|: {max_moment:.1f} in-kips\n'
            'Compression side shown\n'
            'Red marker = Plastic Hinge',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    return fig, ax

def plot_axial_diagram(CORD, ECON, NCT, NE, forces, load_factor):
    """Plot axial force diagram"""
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Find plot limits
    x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
    y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    margin = max(x_range, y_range) * 0.15
    
    # Find max force for scaling
    max_force = np.max(np.abs(forces))
    if max_force > 0:
        force_scale = margin * 0.5 / max_force
    else:
        force_scale = 1.0
    
    # Draw elements and axial force diagrams
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        x1, y1 = CORD[n1, 0], CORD[n1, 1]
        x2, y2 = CORD[n2, 0], CORD[n2, 1]
        
        # Draw element
        ax.plot([x1, x2], [y1, y2], 'b-', linewidth=2.5, zorder=1)
        
        N = forces[el]  # Axial force (constant along element)
        
        # Element direction
        dx = x2 - x1
        dy = y2 - y1
        L = np.sqrt(dx**2 + dy**2)
        
        # Normal direction
        nx = -dy / L
        ny = dx / L
        
        # Draw constant axial force diagram
        offset = N * force_scale
        x_diag = [x1 + offset * nx, x2 + offset * nx]
        y_diag = [y1 + offset * ny, y2 + offset * ny]

        x_full = [x1, x_diag[0], x_diag[1], x2]
        y_full = [y1, y_diag[0], y_diag[1], y2]

        color = 'lightcoral' if N > 0 else 'lightblue'
        ax.fill(x_full, y_full, color=color, alpha=0.5, edgecolor='green', linewidth=1.5)

        # Force value at midpoint offset
        mid_x = (x1 + x2) / 2 + offset * nx * 1.3
        mid_y = (y1 + y2) / 2 + offset * ny * 1.3
        label = f'{N:.1f}\n({"T" if N > 0 else "C"})'
        ax.text(mid_x, mid_y, label, fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # Draw nodes only (no node-number labels on force diagrams)
    for nd in range(NCT):
        ax.plot(CORD[nd, 0], CORD[nd, 1], 'ko', markersize=8, zorder=3)

    ax.set_xlim(x_min - margin * 1.5, x_max + margin * 1.5)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin * 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(r'$X$ (in)', fontsize=12)
    ax.set_ylabel(r'$Y$ (in)', fontsize=12)
    ax.set_title(rf'Axial Force Diagram ($\lambda$ = {load_factor:.3f})',
                 fontsize=14, fontweight='bold')

    ax.text(0.02, 0.98,
            f'Peak |N|: {max_force:.1f} kips\n'
            'Red = Tension (T)\nBlue = Compression (C)',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    return fig, ax

def plot_shear_diagram(CORD, ECON, NCT, NE, moments, load_factor):
    """Plot shear force diagram (derived from moments)"""
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Find plot limits
    x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
    y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    margin = max(x_range, y_range) * 0.15
    
    # V(x) = -dM/dx = -(Mi + Mj)/L  (constant along element)
    shears = np.zeros(NE)
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        dx = CORD[n2, 0] - CORD[n1, 0]
        dy = CORD[n2, 1] - CORD[n1, 1]
        L  = np.sqrt(dx**2 + dy**2)
        shears[el] = -(moments[el, 0] + moments[el, 1]) / L

    max_shear = np.max(np.abs(shears))
    shear_scale = margin * 0.5 / max_shear if max_shear > 0 else 1.0

    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        x1, y1 = CORD[n1, 0], CORD[n1, 1]
        x2, y2 = CORD[n2, 0], CORD[n2, 1]

        ax.plot([x1, x2], [y1, y2], 'b-', linewidth=2.5, zorder=1)

        V  = shears[el]
        dx = x2 - x1;  dy = y2 - y1
        L  = np.sqrt(dx**2 + dy**2)
        nx = -dy / L;  ny = dx / L

        offset = V * shear_scale
        x_full = [x1, x1 + offset*nx, x2 + offset*nx, x2]
        y_full = [y1, y1 + offset*ny, y2 + offset*ny, y2]
        ax.fill(x_full, y_full, color='lightgray', alpha=0.6,
                edgecolor='darkgreen', linewidth=1.5)

        mid_x = (x1 + x2) / 2 + offset * nx * 1.3
        mid_y = (y1 + y2) / 2 + offset * ny * 1.3
        ax.text(mid_x, mid_y, f'{V:.1f}', fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # Draw nodes only (no node-number labels on force diagrams)
    for nd in range(NCT):
        ax.plot(CORD[nd, 0], CORD[nd, 1], 'ko', markersize=8, zorder=3)

    ax.set_xlim(x_min - margin * 1.5, x_max + margin * 1.5)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin * 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel(r'$X$ (in)', fontsize=12)
    ax.set_ylabel(r'$Y$ (in)', fontsize=12)
    ax.set_title(rf'Shear Force Diagram ($\lambda$ = {load_factor:.3f})',
                 fontsize=14, fontweight='bold')

    ax.text(0.02, 0.98,
            f'Peak |V|: {max_shear:.1f} kips\n'
            r'$V = -\mathrm{d}M/\mathrm{d}x = -(M_i+M_j)/L$',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    return fig, ax

def plot_load_displacement(disp_history, load_factors, hinges, CORD, NCT, DI, DLMT, title):
    """
    Plot load factor λ vs. maximum cumulative nodal displacement.

    Each hinge formation is marked.  If DI is available a vertical dashed line
    marks the hinge count at which the degree of static indeterminacy is
    exhausted (the theoretical mechanism threshold).  If DLMT is available a
    horizontal dashed line marks the displacement limit.
    """
    if len(disp_history) == 0:
        return None, None

    fig, ax = plt.subplots(figsize=(10, 7))

    # Maximum resultant nodal displacement at each stage
    max_disp_vals = []
    for disp in disp_history:
        # Resultant translation magnitude at each node
        results = np.sqrt(disp[:, 0]**2 + disp[:, 1]**2)
        max_disp_vals.append(np.max(results))

    lf_arr  = np.array(load_factors)
    md_arr  = np.array(max_disp_vals)

    # Main curve — add origin
    lf_plot  = np.concatenate([[0.0], lf_arr])
    md_plot  = np.concatenate([[0.0], md_arr])
    ax.plot(md_plot, lf_plot, 'b-o', linewidth=2, markersize=6,
            label='Load-displacement path')

    # Mark each hinge formation
    mechanism_marked = False
    for idx, (h_num, el_num, nd_num) in enumerate(hinges):
        lf  = lf_arr[idx]
        md  = md_arr[idx]
        is_mechanism = (DI is not None) and (h_num > DI)
        color = 'orange' if is_mechanism else 'red'
        label = None
        if not is_mechanism and idx == 0:
            label = 'Hinge formation'
        if is_mechanism and not mechanism_marked:
            label = f'Post-mechanism hinge (NCYCL > DI={DI})'
            mechanism_marked = True
        ax.plot(md, lf, 's', color=color, markersize=10, zorder=5, label=label)
        ax.annotate(f'H{h_num}', xy=(md, lf),
                    xytext=(md + max(md_arr)*0.02, lf),
                    fontsize=8, color=color, va='center')

    # Displacement limit line — omitted so the pre-collapse detail is clear
    # (DLMT is still reported in the output file header)

    # DI threshold: horizontal dotted line at λ when DI is exhausted
    if DI is not None and len(hinges) > DI:
        lf_di = lf_arr[DI - 1]   # λ at the DI-th hinge (0-indexed: DI-1)
        ax.axhline(y=lf_di, color='purple', linestyle=':', linewidth=1.5,
                   label=rf'$\lambda$ at DI exhaustion ({DI} hinges) = {lf_di:.2f}')

    ax.set_xlabel(r'Max nodal displacement $\|\Delta\|$ (in)', fontsize=12)
    ax.set_ylabel(r'Load factor $\lambda$', fontsize=12)
    ax.set_title(f'{title}\nLoad-Displacement Curve — Progressive Collapse',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0, right=max(md_arr) * 1.05)
    ax.set_ylim(bottom=0)

    return fig, ax


def visualize_frame(output_file):
    """Create all visualizations for frame analysis from output file only"""

    path = './plots/'
    base = os.path.splitext(os.path.basename(output_file))[0]

    if not os.path.exists(path):
        os.mkdir(path)
        print("sub-directory %s created!" % path)
    else:
        print("sub-directory %s already exists" % path)

    # Read all data from output file
    (FN, NCT, NE, E_mod, CORD, DOF, RTYPE, ECON, PM, SMA,
     hinges, disp_history, moment_history, force_history,
     load_factors, liftoff_history, title, DI, DLMT) = read_output_file(output_file)

    plot_title = title if title else base
    print(f"{base}: {NCT} nodes, {NE} elements")
    print(f"Found {len(hinges)} plastic hinges")
    for h_num, el_num, nd_num in hinges:
        print(f"  Hinge {h_num:2d}: Element {el_num+1}, Node {nd_num+1}")

    # 1. Original geometry
    fig1, ax1 = plot_frame_geometry(CORD, DOF, RTYPE, ECON, NCT, NE,
                                    title=f"{base} — Geometry")
    plt.tight_layout()
    fname = f'{path}{base}-geometry.pdf'
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    print(f"Saved: {fname}")
    plt.close(fig1)

    # 2. Deformed shapes — one per hinge stage
    for idx, (hinge_info, disp, moments, lf) in enumerate(
            zip(hinges, disp_history, moment_history, load_factors)):
        hinges_so_far = hinges[:idx+1]
        liftoff = liftoff_history[idx] if idx < len(liftoff_history) else set()
        fig2, ax2 = plot_deformed_shape(CORD, DOF, RTYPE, ECON, NCT, NE,
                                        disp, moments, lf, SMA, E_mod,
                                        scale=None, hinge_info=hinges_so_far,
                                        liftoff_nodes=liftoff)
        plt.tight_layout()
        fname = f'{path}{base}-deformed_hinge_{idx+1:02d}.pdf'
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"Saved: {fname}")
        plt.close(fig2)

    # 3. Moment diagrams — one per hinge stage
    for idx, (hinge_info, moments, lf) in enumerate(zip(hinges, moment_history, load_factors)):
        hinges_so_far = hinges[:idx+1]
        fig3, ax3 = plot_moment_diagram(CORD, ECON, NCT, NE, moments, PM, lf,
                                        hinges=hinges_so_far)
        plt.tight_layout()
        fname = f'{path}{base}-moments_hinge_{idx+1:02d}.pdf'
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"Saved: {fname}")
        plt.close(fig3)

    # 4. Axial force diagram — final state
    n_hinges = len(hinges)
    if len(force_history) > 0:
        fig4, ax4 = plot_axial_diagram(CORD, ECON, NCT, NE,
                                       force_history[-1], load_factors[-1])
        plt.tight_layout()
        fname = f'{path}{base}-axial_hinge_{n_hinges:02d}.pdf'
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"Saved: {fname}")
        plt.close(fig4)

    # 5. Shear force diagram — final state
    if len(moment_history) > 0:
        fig5, ax5 = plot_shear_diagram(CORD, ECON, NCT, NE,
                                       moment_history[-1], load_factors[-1])
        plt.tight_layout()
        fname = f'{path}{base}-shear_hinge_{n_hinges:02d}.pdf'
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"Saved: {fname}")
        plt.close(fig5)

    # 6. Summary: 2×2 progressive collapse panels
    if len(hinges) > 0:
        n_stages = min(4, len(hinges))
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()

        for idx in range(n_stages):
            ax = axes[idx]

            disp    = disp_history[idx]
            moments = moment_history[idx]
            lf      = load_factors[idx]
            hinges_so_far = hinges[:idx+1]

            x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
            y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
            x_range, y_range = x_max - x_min, y_max - y_min
            margin = max(x_range, y_range) * 0.15
            frame_size = max(x_range, y_range)

            max_disp = np.max(np.abs(disp[:, :2]))
            scale = frame_size * 0.08 / max_disp if max_disp > 0 else 1.0

            # Draw original (dashed)
            for el in range(NE):
                n1, n2 = ECON[el, 0], ECON[el, 1]
                ax.plot([CORD[n1, 0], CORD[n2, 0]], [CORD[n1, 1], CORD[n2, 1]],
                        'b--', linewidth=1, alpha=0.3)

            # Draw curved deformed shape
            curves = compute_element_curves(CORD, ECON, SMA, E_mod,
                                            disp, moments, NE, frame_size, scale=scale)
            for curve in curves:
                ax.plot(curve[:, 0], curve[:, 1], 'r-', linewidth=2)

            # Hinge markers on deformed curves
            for h_num, el_num, nd_num in hinges_so_far:
                curve = curves[el_num]
                x_h = curve[0, 0] if nd_num == ECON[el_num, 0] else curve[-1, 0]
                y_h = curve[0, 1] if nd_num == ECON[el_num, 0] else curve[-1, 1]
                circle = Circle((x_h, y_h), margin * 0.05,
                                facecolor='white', edgecolor='red', linewidth=2)
                ax.add_patch(circle)

            ax.set_xlim(x_min - margin, x_max + margin)
            ax.set_ylim(y_min - margin, y_max + margin)
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.set_title(rf'After Hinge {idx+1:02d} ($\lambda$ = {lf:.3f})',
                         fontsize=11, fontweight='bold')

        plt.suptitle(f'{base} — Progressive Collapse', fontsize=14, fontweight='bold')
        plt.tight_layout()
        fname = f'{path}{base}-summary.pdf'
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        print(f"Saved: {fname}")
        plt.close(fig)

    # 7. Load-displacement curve
    if len(hinges) > 0:
        fig7, ax7 = plot_load_displacement(
            disp_history, load_factors, hinges, CORD, NCT, DI, DLMT, plot_title)
        if fig7 is not None:
            plt.tight_layout()
            fname = f'{path}{base}-load_displacement.pdf'
            plt.savefig(fname, dpi=150, bbox_inches='tight')
            print(f"Saved: {fname}")
            plt.close(fig7)

    print("\nVisualization complete!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python epframe_viz.py output_file")
        print("\nThe output file contains all necessary data (geometry + results)")
        sys.exit(1)
    
    visualize_frame(sys.argv[1])
