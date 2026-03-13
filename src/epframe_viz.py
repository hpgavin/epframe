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
from matplotlib.patches import Circle, Polygon, FancyBboxPatch
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

    FN = 0
    NCT = 0
    NE = 0
    E_mod = 0.0
    CORD = None
    ECON = None
    PM = None
    DOF = None
    RTYPE = None

    # --- Pass 1: parse header / input echo ---
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if 'FRAME NO' in line:
            parts = line.split()
            for j, p in enumerate(parts):
                if p == 'NO':
                    FN = int(parts[j+1])
                    break

        if 'NUMBER OF NODES'    in line: NCT   = int(line.split()[-1])
        if 'NUMBER OF ELEMENTS' in line: NE    = int(line.split()[-1])
        if 'MOD OF ELASTICITY'  in line: E_mod = float(line.split()[-1])

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

        if 'DATA FOR ELEMENTS' in line:
            i += 2
            ECON = np.zeros((NE, 2), dtype=int)
            PM   = np.zeros(NE)

            for el in range(NE):
                i += 1
                data_line = lines[i].replace('%', '').strip()
                parts = data_line.split()
                elem_num = int(parts[0]) - 1
                ECON[elem_num, 0] = int(parts[1]) - 1
                ECON[elem_num, 1] = int(parts[2]) - 1
                PM[elem_num]      = float(parts[5])

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
            while i < len(lines) and 'TENSION' not in lines[i]:
                i += 1
            i += 1

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

    return (FN, NCT, NE, E_mod, CORD, DOF, RTYPE, ECON, PM,
            hinges, disp_history, moment_history, force_history,
            load_factors, liftoff_history)

def draw_support(ax, x, y, DOF, RTYPE, size, lifted=False):
    """
    Draw support symbol based on DOF flags and reaction types.

    DOF[i]   : 0 = constrained, 1 = free
    RTYPE[i] : 0 = free, 1 = bidirectional, 2 = positive-only, 3 = negative-only
    lifted   : True when a unidirectional support has lifted off (draw gap symbol)

    Support classification:
      Fully fixed       : RTYPE == [1,1,1]  (no free DOFs)
      Pinned            : RTYPE[0]==1, RTYPE[1]==1, RTYPE[2]==0
      Bidirectional roller: one translation constrained bidirectionally, rotation free
      Unidirectional roller: one translation constrained POS or NEG, rotation free
    """
    rx, ry, rz = RTYPE[0], RTYPE[1], RTYPE[2]

    if lifted:
        # Lifted-off unidirectional support: dashed outline + gap marker
        circle = plt.Circle((x, y - size*0.5), size*0.35,
                             facecolor='none', edgecolor='gray',
                             linewidth=1.5, linestyle='--', zorder=2)
        ax.add_patch(circle)
        ax.plot([x - size, x + size], [y - size*0.9, y - size*0.9],
                'k--', linewidth=1.5, alpha=0.5)
        # Gap lines
        for dx_g in [-size*0.4, 0, size*0.4]:
            ax.plot([x + dx_g, x + dx_g - size*0.15],
                    [y - size*0.9, y - size*1.2], 'k-', linewidth=1, alpha=0.4)
        return

    if rx == 1 and ry == 1 and rz == 1:
        # Fully fixed — filled triangle + hatching
        tri = Polygon([(x, y), (x - size, y - size*1.2), (x + size, y - size*1.2)],
                      closed=True, facecolor='gray', edgecolor='black', linewidth=1.5)
        ax.add_patch(tri)
        for j in range(5):
            xh = x - size + j * size * 0.5
            ax.plot([xh, xh - size*0.3], [y - size*1.2, y - size*1.6], 'k-', linewidth=1)

    elif rx == 1 and ry == 1 and rz == 0:
        # Pinned — hollow triangle + ground line
        tri = Polygon([(x, y), (x - size, y - size*1.2), (x + size, y - size*1.2)],
                      closed=True, facecolor='white', edgecolor='black', linewidth=1.5)
        ax.add_patch(tri)
        ax.plot([x - size*1.2, x + size*1.2], [y - size*1.2, y - size*1.2], 'k-', linewidth=2)

    elif ry == 1 and rx != 1:
        # Standard bidirectional Y roller — circle + horizontal line
        circle = plt.Circle((x, y - size*0.4), size*0.35,
                             facecolor='white', edgecolor='black', linewidth=1.5)
        ax.add_patch(circle)
        ax.plot([x - size, x + size], [y - size*0.8, y - size*0.8], 'k-', linewidth=2)

    elif ry == 2:
        # Unidirectional Y+ roller (reaction force upward, resists downward displacement)
        # Symbol: roller circle with upward arrow
        circle = plt.Circle((x, y - size*0.4), size*0.35,
                             facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(circle)
        ax.plot([x - size, x + size], [y - size*0.8, y - size*0.8], 'k-', linewidth=2)
        ax.annotate('', xy=(x, y + size*0.1), xytext=(x, y - size*0.8),
                    arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))

    elif ry == 3:
        # Unidirectional Y− roller (reaction force downward, resists upward displacement)
        # Symbol: roller circle with downward arrow
        circle = plt.Circle((x, y - size*0.4), size*0.35,
                             facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(circle)
        ax.plot([x - size, x + size], [y - size*0.8, y - size*0.8], 'k-', linewidth=2)
        ax.annotate('', xy=(x, y - size*0.8), xytext=(x, y + size*0.1),
                    arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))

    elif rx == 2:
        # Unidirectional X+ roller (reaction force rightward)
        circle = plt.Circle((x - size*0.4, y), size*0.35,
                             facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(circle)
        ax.plot([x - size*0.8, x - size*0.8], [y - size, y + size], 'k-', linewidth=2)
        ax.annotate('', xy=(x + size*0.1, y), xytext=(x - size*0.8, y),
                    arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))

    elif rx == 3:
        # Unidirectional X− roller (reaction force leftward)
        circle = plt.Circle((x + size*0.4, y), size*0.35,
                             facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(circle)
        ax.plot([x + size*0.8, x + size*0.8], [y - size, y + size], 'k-', linewidth=2)
        ax.annotate('', xy=(x - size*0.1, y), xytext=(x + size*0.8, y),
                    arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))

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
    
    # Draw nodes and supports
    for nd in range(NCT):
        x, y = CORD[nd, 0], CORD[nd, 1]
        
        # Draw support if any DOF is constrained
        if np.any(RTYPE[nd] > 0):
            draw_support(ax, x, y, DOF[nd], RTYPE[nd], margin * 0.15)
        
        # Draw node
        ax.plot(x, y, 'ko', markersize=8, zorder=3)
        ax.text(x, y + margin * 0.12, f'{nd+1}', fontsize=11, fontweight='bold',
                ha='center', va='bottom', zorder=4)
    
    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (in)', fontsize=12)
    ax.set_ylabel('Y (in)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    return fig, ax

def plot_deformed_shape(CORD, DOF, RTYPE, ECON, NCT, NE, disp, load_factor,
                        scale=None, hinge_info=None, liftoff_nodes=None):
    """Plot deformed shape with optional hinge markers"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Find plot limits
    x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
    y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
    x_range = x_max - x_min
    y_range = y_max - y_min
    margin = max(x_range, y_range) * 0.15
    
    # Auto-calculate scale if not provided
    if scale is None:
        max_disp = np.max(np.abs(disp[:, :2]))
        if max_disp > 0:
            scale = max(x_range, y_range) * 0.08 / max_disp
        else:
            scale = 1.0
    
    # Draw original shape (dashed)
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        ax.plot([CORD[n1, 0], CORD[n2, 0]], [CORD[n1, 1], CORD[n2, 1]],
                'b--', linewidth=1.5, alpha=0.4, label='Original' if el == 0 else '')
    
    # Draw deformed shape
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        x1_def = CORD[n1, 0] + scale * disp[n1, 0]
        y1_def = CORD[n1, 1] + scale * disp[n1, 1]
        x2_def = CORD[n2, 0] + scale * disp[n2, 0]
        y2_def = CORD[n2, 1] + scale * disp[n2, 1]
        ax.plot([x1_def, x2_def], [y1_def, y2_def], 'r-', linewidth=2.5,
                label='Deformed' if el == 0 else '')
    
    # Draw nodes
    for nd in range(NCT):
        x_def = CORD[nd, 0] + scale * disp[nd, 0]
        y_def = CORD[nd, 1] + scale * disp[nd, 1]
        ax.plot(x_def, y_def, 'ro', markersize=6, zorder=3)
        ax.text(x_def, y_def + margin * 0.08, f'{nd+1}', fontsize=10,
                ha='center', va='bottom', zorder=4)
    
    # Draw supports at original positions; lifted-off supports shown with gap symbol
    _coord_name = {0: 'X', 1: 'Y', 2: 'Z'}
    for nd in range(NCT):
        if np.any(RTYPE[nd] > 0):
            # A support is lifted if any of its unidirectional DOFs are in liftoff_nodes
            lifted = liftoff_nodes is not None and any(
                (nd, _coord_name[c]) in liftoff_nodes
                for c in range(3) if RTYPE[nd, c] in (2, 3)
            )
            draw_support(ax, CORD[nd, 0], CORD[nd, 1], DOF[nd], RTYPE[nd],
                         margin * 0.12, lifted=lifted)
    
    # Draw plastic hinges
    if hinge_info:
        for h_num, el_num, nd_num in hinge_info:
            x_h = CORD[nd_num, 0] + scale * disp[nd_num, 0]
            y_h = CORD[nd_num, 1] + scale * disp[nd_num, 1]
            circle = Circle((x_h, y_h), margin * 0.06,
                           facecolor='white', edgecolor='red', linewidth=2.5, zorder=10)
            ax.add_patch(circle)
            ax.text(x_h, y_h, f'{h_num}', fontsize=9, fontweight='bold',
                   color='red', ha='center', va='center', zorder=11)
    
    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (in)', fontsize=12)
    ax.set_ylabel('Y (in)', fontsize=12)
    ax.set_title(f'Deformed Shape (Load Factor λ = {load_factor:.3f}, Scale = {scale:.1f}x)',
                fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    
    # Add displacement info
    max_x = np.max(np.abs(disp[:, 0]))
    max_y = np.max(np.abs(disp[:, 1]))
    max_r = np.max(np.abs(disp[:, 2]))
    info_text = f'Max displacements:\nX: {max_x:.4f} in\nY: {max_y:.4f} in\nRot: {max_r:.5f} rad'
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
        
        M1 = moments[el, 0]  # Near-end moment
        M2 = moments[el, 1]  # Far-end moment
        
        # Element direction
        dx = x2 - x1
        dy = y2 - y1
        L = np.sqrt(dx**2 + dy**2)
        
        # Normal direction (perpendicular)
        nx = -dy / L
        ny = dx / L
        
        # Create moment diagram (linear variation)
        n_points = 20
        x_diag = []
        y_diag = []
        
        for j in range(n_points + 1):
            s = j / n_points
            x_pt = x1 + s * dx
            y_pt = y1 + s * dy
            M_pt = M1 * (1 - s) + M2 * s
            
            x_diag.append(x_pt + M_pt * moment_scale * nx)
            y_diag.append(y_pt + M_pt * moment_scale * ny)
        
        # Draw filled moment diagram - use grey for all
        x_full = [x1] + x_diag + [x2]
        y_full = [y1] + y_diag + [y2]
        
        ax.fill(x_full, y_full, color='lightgray', alpha=0.6, edgecolor='blue', linewidth=1.5)
        
        # Add moment values at ends
        ax.text(x1 + M1*moment_scale*nx*1.2, y1 + M1*moment_scale*ny*1.2,
                f'{M1:.0f}', fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(x2 + M2*moment_scale*nx*1.2, y2 + M2*moment_scale*ny*1.2,
                f'{M2:.0f}', fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # Label element
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        ax.text(mid_x, mid_y, f'E{el+1}', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                ha='center', va='center', zorder=5)
        
        # Mark plastic hinges
        if hinges:
            for h_num, el_num, nd_num in hinges:
                if el_num == el:
                    if nd_num == n1:
                        x_h, y_h = x1, y1
                    else:
                        x_h, y_h = x2, y2
                    
                    circle = Circle((x_h, y_h), margin * 0.06,
                                   facecolor='white', edgecolor='red', linewidth=2.5, zorder=10)
                    ax.add_patch(circle)
                    ax.text(x_h, y_h - margin * 0.12, f'H{h_num}',
                            fontsize=9, fontweight='bold', color='red',
                            ha='center', va='top')
    
    # Draw nodes
    for nd in range(NCT):
        ax.plot(CORD[nd, 0], CORD[nd, 1], 'ko', markersize=8, zorder=3)
        ax.text(CORD[nd, 0], CORD[nd, 1] + margin * 0.08, f'{nd+1}',
                fontsize=11, fontweight='bold', ha='center', va='bottom', zorder=4)
    
    ax.set_xlim(x_min - margin * 1.5, x_max + margin * 1.5)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin * 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (in)', fontsize=12)
    ax.set_ylabel('Y (in)', fontsize=12)
    ax.set_title(f'Bending Moment Diagram (Load Factor λ = {load_factor:.3f})',
                fontsize=14, fontweight='bold')
    
    # Add legend
    ax.text(0.02, 0.98,
            f'Max Moment: {max_moment:.1f} in-kips\n' +
            f'Plastic Moment: {PM[0]:.1f} in-kips\n' +
            f'Tension side shown\nRed marker = Plastic Hinge',
            transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
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
        
        # Color: red for tension, blue for compression
        color = 'lightcoral' if N > 0 else 'lightblue'
        ax.fill(x_full, y_full, color=color, alpha=0.5, edgecolor='green', linewidth=1.5)
        
        # Add force value at midpoint
        mid_x = (x1 + x2) / 2 + offset * nx * 1.3
        mid_y = (y1 + y2) / 2 + offset * ny * 1.3
        label = f'{N:.1f}\n({"T" if N > 0 else "C"})'
        ax.text(mid_x, mid_y, label, fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # Label element
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        ax.text(mid_x, mid_y, f'E{el+1}', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                ha='center', va='center', zorder=5)
    
    # Draw nodes
    for nd in range(NCT):
        ax.plot(CORD[nd, 0], CORD[nd, 1], 'ko', markersize=8, zorder=3)
        ax.text(CORD[nd, 0], CORD[nd, 1] + margin * 0.08, f'{nd+1}',
                fontsize=11, fontweight='bold', ha='center', va='bottom', zorder=4)
    
    ax.set_xlim(x_min - margin * 1.5, x_max + margin * 1.5)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin * 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (in)', fontsize=12)
    ax.set_ylabel('Y (in)', fontsize=12)
    ax.set_title(f'Axial Force Diagram (Load Factor λ = {load_factor:.3f})',
                fontsize=14, fontweight='bold')
    
    # Add legend
    ax.text(0.02, 0.98,
            f'Max Axial: {max_force:.1f} kips\n' +
            f'Red = Tension (T)\nBlue = Compression (C)',
            transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
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
    
    # Calculate shear forces from moment equilibrium: V = (M1 + M2) / L
    shears = np.zeros(NE)
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        dx = CORD[n2, 0] - CORD[n1, 0]
        dy = CORD[n2, 1] - CORD[n1, 1]
        L = np.sqrt(dx**2 + dy**2)
        shears[el] = (moments[el, 0] + moments[el, 1]) / L
    
    # Find max shear for scaling
    max_shear = np.max(np.abs(shears))
    if max_shear > 0:
        shear_scale = margin * 0.5 / max_shear
    else:
        shear_scale = 1.0
    
    # Draw elements and shear diagrams
    for el in range(NE):
        n1, n2 = ECON[el, 0], ECON[el, 1]
        x1, y1 = CORD[n1, 0], CORD[n1, 1]
        x2, y2 = CORD[n2, 0], CORD[n2, 1]
        
        # Draw element
        ax.plot([x1, x2], [y1, y2], 'b-', linewidth=2.5, zorder=1)
        
        V = shears[el]
        
        # Element direction
        dx = x2 - x1
        dy = y2 - y1
        L = np.sqrt(dx**2 + dy**2)
        
        # Normal direction
        nx = -dy / L
        ny = dx / L
        
        # Draw constant shear diagram
        offset = V * shear_scale
        x_diag = [x1 + offset * nx, x2 + offset * nx]
        y_diag = [y1 + offset * ny, y2 + offset * ny]
        
        x_full = [x1, x_diag[0], x_diag[1], x2]
        y_full = [y1, y_diag[0], y_diag[1], y2]
        
        # Use grey shading for shear
        ax.fill(x_full, y_full, color='lightgray', alpha=0.6, edgecolor='darkgreen', linewidth=1.5)
        
        # Add shear value at midpoint
        mid_x = (x1 + x2) / 2 + offset * nx * 1.3
        mid_y = (y1 + y2) / 2 + offset * ny * 1.3
        ax.text(mid_x, mid_y, f'{V:.1f}', fontsize=9, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        
        # Label element
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        ax.text(mid_x, mid_y, f'E{el+1}', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                ha='center', va='center', zorder=5)
    
    # Draw nodes
    for nd in range(NCT):
        ax.plot(CORD[nd, 0], CORD[nd, 1], 'ko', markersize=8, zorder=3)
        ax.text(CORD[nd, 0], CORD[nd, 1] + margin * 0.08, f'{nd+1}',
                fontsize=11, fontweight='bold', ha='center', va='bottom', zorder=4)
    
    ax.set_xlim(x_min - margin * 1.5, x_max + margin * 1.5)
    ax.set_ylim(y_min - margin * 1.5, y_max + margin * 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('X (in)', fontsize=12)
    ax.set_ylabel('Y (in)', fontsize=12)
    ax.set_title(f'Shear Force Diagram (Load Factor λ = {load_factor:.3f})',
                fontsize=14, fontweight='bold')
    
    # Add legend
    ax.text(0.02, 0.98,
            f'Max Shear: {max_shear:.1f} kips\n' +
            f'V = (M₁ + M₂) / L',
            transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    return fig, ax

def visualize_frame(output_file):
    """Create all visualizations for frame analysis from output file only"""


    path = './plots/'

    # check whether './plots/' already exists
    if not os.path.exists(path):
        os.mkdir(path)
        print("sub-directory %s created!" % path)
    else:
        print("sub-directory %s already exists" % path)
    
    # Read all data from output file
    (FN, NCT, NE, E_mod, CORD, DOF, RTYPE, ECON, PM,
     hinges, disp_history, moment_history, force_history,
     load_factors, liftoff_history) = read_output_file(output_file)
    
    print(f"Frame {FN}: {NCT} nodes, {NE} elements")
    print(f"Found {len(hinges)} plastic hinges")
    for h_num, el_num, nd_num in hinges:
        print(f"  Hinge {h_num}: Element {el_num+1}, Node {nd_num+1}")
    
    # 1. Plot original geometry
    fig1, ax1 = plot_frame_geometry(CORD, DOF, RTYPE, ECON, NCT, NE,
                                    title=f"Frame {FN} - Original Geometry")
    plt.tight_layout()
    plt.savefig(f'{path}frame_{FN}_geometry.pdf', dpi=150, bbox_inches='tight')
    print(f"Saved: {path}frame_{FN}_geometry.pdf")
    plt.close(fig1)
    
    # 2. Plot deformed shapes for each hinge formation
    for idx, (hinge_info, disp, lf) in enumerate(zip(hinges, disp_history, load_factors)):
        hinges_so_far = hinges[:idx+1]
        liftoff = liftoff_history[idx] if idx < len(liftoff_history) else set()
        fig2, ax2 = plot_deformed_shape(CORD, DOF, RTYPE, ECON, NCT, NE, disp, lf,
                                        scale=None, hinge_info=hinges_so_far,
                                        liftoff_nodes=liftoff)
        plt.tight_layout()
        plt.savefig(f'{path}frame_{FN}_deformed_hinge_{idx+1}.pdf', dpi=150, bbox_inches='tight')
        print(f"Saved: {path}frame_{FN}_deformed_hinge_{idx+1}.pdf")
        plt.close(fig2)
    
    # 3. Plot moment diagrams for each hinge formation
    for idx, (hinge_info, moments, lf) in enumerate(zip(hinges, moment_history, load_factors)):
        hinges_so_far = hinges[:idx+1]
        fig3, ax3 = plot_moment_diagram(CORD, ECON, NCT, NE, moments, PM, lf,
                                        hinges=hinges_so_far)
        plt.tight_layout()
        plt.savefig(f'{path}frame_{FN}_moments_hinge_{idx+1}.pdf', dpi=150, bbox_inches='tight')
        print(f"Saved: {path}frame_{FN}_moments_hinge_{idx+1}.pdf")
        plt.close(fig3)
    
    # 4. Plot axial force diagrams for final state
    if len(force_history) > 0:
        fig4, ax4 = plot_axial_diagram(CORD, ECON, NCT, NE, force_history[-1], load_factors[-1])
        plt.tight_layout()
        plt.savefig(f'{path}frame_{FN}_axial_{idx+1}.pdf', dpi=150, bbox_inches='tight')
        print(f"Saved: {path}frame_{FN}_axial_{idx+1}.pdf")
        plt.close(fig4)
    
    # 5. Plot shear force diagrams for final state
    if len(moment_history) > 0:
        fig5, ax5 = plot_shear_diagram(CORD, ECON, NCT, NE, moment_history[-1], load_factors[-1])
        plt.tight_layout()
        plt.savefig(f'{path}frame_{FN}_shear_{idx+1}.pdf', dpi=150, bbox_inches='tight')
        print(f"Saved: {path}frame_{FN}_shear_{idx+1}.pdf")
        plt.close(fig5)
    
    # 6. Create summary plot showing progressive collapse
    if len(hinges) > 0:
        n_stages = min(4, len(hinges))
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        for idx in range(n_stages):
            ax = axes[idx]
            
            disp = disp_history[idx]
            lf = load_factors[idx]
            hinges_so_far = hinges[:idx+1]
            
            # Find plot limits
            x_min, x_max = CORD[:, 0].min(), CORD[:, 0].max()
            y_min, y_max = CORD[:, 1].min(), CORD[:, 1].max()
            x_range, y_range = x_max - x_min, y_max - y_min
            margin = max(x_range, y_range) * 0.15
            
            # Auto-scale
            max_disp = np.max(np.abs(disp[:, :2]))
            if max_disp > 0:
                scale = max(x_range, y_range) * 0.08 / max_disp
            else:
                scale = 1.0
            
            # Draw original (dashed)
            for el in range(NE):
                n1, n2 = ECON[el, 0], ECON[el, 1]
                ax.plot([CORD[n1, 0], CORD[n2, 0]], [CORD[n1, 1], CORD[n2, 1]],
                        'b--', linewidth=1, alpha=0.3)
            
            # Draw deformed (solid)
            for el in range(NE):
                n1, n2 = ECON[el, 0], ECON[el, 1]
                x1_def = CORD[n1, 0] + scale * disp[n1, 0]
                y1_def = CORD[n1, 1] + scale * disp[n1, 1]
                x2_def = CORD[n2, 0] + scale * disp[n2, 0]
                y2_def = CORD[n2, 1] + scale * disp[n2, 1]
                ax.plot([x1_def, x2_def], [y1_def, y2_def], 'r-', linewidth=2)
            
            # Draw hinges
            for h_num, el_num, nd_num in hinges_so_far:
                x_h = CORD[nd_num, 0] + scale * disp[nd_num, 0]
                y_h = CORD[nd_num, 1] + scale * disp[nd_num, 1]
                circle = Circle((x_h, y_h), margin * 0.05,
                               facecolor='white', edgecolor='red', linewidth=2)
                ax.add_patch(circle)
            
            ax.set_xlim(x_min - margin, x_max + margin)
            ax.set_ylim(y_min - margin, y_max + margin)
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.set_title(f'After Hinge {idx+1} (λ = {lf:.3f})', fontsize=11, fontweight='bold')
        
        plt.suptitle(f'Frame {FN} - Progressive Collapse Analysis',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{path}frame_{FN}_summary.pdf', dpi=150, bbox_inches='tight')
        print(f"Saved: {path}frame_{FN}_summary.pdf")
        plt.close(fig)
    
    print("\nVisualization complete!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python epframe_viz.py output_file")
        print("\nThe output file contains all necessary data (geometry + results)")
        sys.exit(1)
    
    visualize_frame(sys.argv[1])
