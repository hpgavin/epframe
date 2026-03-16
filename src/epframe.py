#!/usr/bin/env python3
"""
ELASTIC PLASTIC ANALYSIS OF A PLANE FRAME WITH ONE-WAY REACTIONS
Based on epframe.py with QP solver for inequality constraints

One-way reaction types:
  '0' or 0   : No reaction (free DOF)
  '*' or 1   : Bi-directional reaction (standard support)
  '+' or 2   : Reaction only in positive direction
  '-' or 3   : Reaction only in negative direction
  
Nomenclature:
    Node (N) - connection point in the structure (0-indexed internally, 1-indexed in I/O)
    Element (E) - structural member connecting two nodes (0-indexed internally, 1-indexed in I/O)
"""

import numpy as np
import sys
import warnings
from math import sqrt
from scipy.linalg import solve
from scipy.linalg import LinAlgWarning
from numpy.linalg import LinAlgError

def get_csv_header(NCT, NE):
    """Generate CSV header row"""
    headers = ['NCYCL', 'EL', 'NH', 'CLF']
    # CD - cumulative displacement - 3 per node (X, Y, Rotation)
    for nd in range(NCT):
        headers.append(f'CD{nd+1}_X')
        headers.append(f'CD{nd+1}_Y')
        headers.append(f'CD{nd+1}_R')
    # CM - cumulative moment - 2 per element (moment at each end)
    for el in range(NE):
        headers.append(f'CM{el+1}_1')
        headers.append(f'CM{el+1}_2')
    # CT - cumulative tension - 1 per element (axial force)
    for el in range(NE):
        headers.append(f'CT{el+1}')
    # Active status - 3 per node
    for nd in range(NCT):
        headers.append(f'ACT{nd+1}_X')
        headers.append(f'ACT{nd+1}_Y')
        headers.append(f'ACT{nd+1}_R')
    return ','.join(headers)

def write_csv_row(csv_fp, NCYCL, EL, NH, CLF, CD, CM, CT, DOF, RTYPE, active, NCT, NE):
    """Write a data row to CSV file.
    Moments and rotations are negated at output to convert from the internal
    CW-positive element convention to the standard CCW-positive convention."""
    values = [f'{NCYCL}', f'{EL}', f'{NH}', f'{CLF:.6E}']

    # CD — negate rotation (dof=2) to convert to CCW-positive
    NN = 0
    for nd in range(NCT):
        for dof in range(3):
            if DOF[nd, dof]:
                val = -CD[NN] if dof == 2 else CD[NN]
                values.append(f'{val:.6E}')
                NN += 1
            else:
                values.append(f'{0.0:.6E}')

    # CM — negate all moments to convert to CCW-positive
    for i in range(2*NE):
        values.append(f'{-CM[i]:.6E}')

    # CT — axial forces unchanged
    for i in range(NE):
        values.append(f'{CT[i]:.6E}')

    # Active status
    NN = 0
    for nd in range(NCT):
        for dof in range(3):
            if DOF[nd, dof]:
                values.append(f'{int(active[NN])}')
                NN += 1
            else:
                values.append('0')

    csv_fp.write(','.join(values) + '\n')

def parse_fixity(fixity_str):
    """
    Parse fixity string to DOF flags and reaction types
    
    Examples:
        'X Y Z' or 'X* Y* Z*' -> DOF=[0,0,0], RTYPE=[1,1,1]
        'X Y *' or 'X* Y* 0' -> DOF=[0,0,1], RTYPE=[1,1,0]
        'X+ Y Z' -> DOF=[0,0,0], RTYPE=[2,1,1]
        '0 0 0' or '* * *' -> DOF=[1,1,1], RTYPE=[0,0,0]
    
    Returns:
        DOF: [DFX, DFY, DFZ] where 1=free, 0=has constraint
        RTYPE: [RTX, RTY, RTZ] where 0=free, 1=bi, 2=pos, 3=neg
    """
    fixity_str = str(fixity_str).upper().strip()
    
    # Initialize as free
    DOF = [1, 1, 1]
    RTYPE = [0, 0, 0]
    
    # Split by whitespace
    parts = fixity_str.split()
    
    for part in parts:
        part = part.strip()
        if not part:
            continue
        
        # Determine which coordinate
        coord_idx = None
        constraint_char = None
        
        if 'X' in part:
            coord_idx = 0
            constraint_char = part.replace('X', '')
        elif 'Y' in part:
            coord_idx = 1
            constraint_char = part.replace('Y', '')
        elif 'Z' in part:
            coord_idx = 2
            constraint_char = part.replace('Z', '')
        else:
            continue
        
        if coord_idx is None:
            continue
        
        # Parse constraint type
        if constraint_char == '' or constraint_char == '*' or constraint_char == '1':
            DOF[coord_idx] = 0
            RTYPE[coord_idx] = 1  # Bidirectional
        elif constraint_char == '+' or constraint_char == '2':
            DOF[coord_idx] = 0
            RTYPE[coord_idx] = 2  # Positive only
        elif constraint_char == '-' or constraint_char == '3':
            DOF[coord_idx] = 0
            RTYPE[coord_idx] = 3  # Negative only
        elif constraint_char == '0':
            DOF[coord_idx] = 1
            RTYPE[coord_idx] = 0  # Free
    
    # Ensure consistency:
    #   Bidirectional (RTYPE==1): permanently constrained, excluded from global DOF vector (DOF=0)
    #   One-way (RTYPE==2 or 3): conditionally constrained; MUST stay in global DOF vector (DOF=1)
    #                             so the active-set solver can enforce or release the constraint
    #   Free (RTYPE==0): unconstrained, included in global DOF vector (DOF=1)
    for i in range(3):
        if RTYPE[i] == 1:
            DOF[i] = 0   # permanently removed from displacement vector
        else:
            DOF[i] = 1   # in displacement vector (free or conditionally constrained)
    
    return DOF, RTYPE

def read_input_file(filename):
    """Read input file with one-way reaction support"""
    with open(filename, 'r') as f:
        lines = []
        for line in f:
            line = line.split('#')[0].strip()
            if line:
                lines.append(line)
    
    idx = 0
    title = lines[idx].strip()   # first line is a free-text title, not a number
    idx += 1
    
    parts = lines[idx].split()
    NCT, NE, E, Fy = int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3])
    idx += 1
    
    # Node data
    CORD = np.zeros((NCT, 2))
    DOF = np.zeros((NCT, 3), dtype=int)
    RTYPE = np.zeros((NCT, 3), dtype=int)
    
    for i in range(NCT):
        parts = lines[idx].split()
        node_num = int(parts[0]) - 1
        CORD[node_num, 0] = float(parts[1])
        CORD[node_num, 1] = float(parts[2])
        
        # Check format - try to distinguish old numeric (0/1) from new string format
        # Old format: three separate integers (0 or 1)
        # New format: strings with X, Y, Z, *, +, -, or 0
        try:
            # Try parsing as old numeric format
            if len(parts) >= 6:
                rx = int(parts[3])
                ry = int(parts[4])
                rz = int(parts[5])
                # Successfully parsed as integers - use old format
                DOF[node_num, 0] = 1 - rx
                DOF[node_num, 1] = 1 - ry
                DOF[node_num, 2] = 1 - rz
                RTYPE[node_num, :] = [1 if d == 0 else 0 for d in DOF[node_num, :]]
            else:
                raise ValueError("Not enough parts for old format")
        except ValueError:
            # New format - parse as fixity string
            fixity_str = ' '.join(parts[3:])
            dof, rtype = parse_fixity(fixity_str)
            DOF[node_num, :] = dof
            RTYPE[node_num, :] = rtype
        
        idx += 1
    
    # Element data
    ECON = np.zeros((NE, 2), dtype=int)
    SMA = np.zeros(NE)
    AREA = np.zeros(NE)
    ZS = np.zeros(NE)    # plastic section modulus Z (in^3)
    
    for i in range(NE):
        parts = lines[idx].split()
        elem_num = int(parts[0]) - 1
        ECON[elem_num, 0] = int(parts[1]) - 1
        ECON[elem_num, 1] = int(parts[2]) - 1
        SMA[elem_num]  = float(parts[3])    # second moment of area, I
        AREA[elem_num] = float(parts[4])    # cross section area, A
        ZS[elem_num]   = float(parts[5])    # plastic section modulus, Z
        idx += 1
    
    # Derived section capacities
    PM = ZS * Fy          # plastic moment  Mp = Z * Fy  (in-kips)
    PY = AREA * Fy        # axial yield force  Py = A * Fy  (kips)
    
    # Load data
    LN = int(lines[idx].split()[0])
    idx += 1
    
    loads = []
    for i in range(LN):
        parts = lines[idx].split()
        NH = int(parts[0]) - 1
        FX = float(parts[1])
        FY = float(parts[2])
        FZ = float(parts[3])
        loads.append((NH, FX, FY, FZ))
        idx += 1
    
    return title, NCT, NE, E, Fy, CORD, DOF, RTYPE, ECON, SMA, AREA, ZS, PM, PY, loads

def build_geometric_stiffness(CORD, ECON, DOF, CT, OLEN_elems, NE, NCT, ND):
    """
    Assemble the global geometric stiffness matrix from current axial forces.

    For a Bernoulli-Euler beam element with axial force N and length L the
    local geometric stiffness matrix is (eq 75 in StructuralElements.pdf):

        K̄_g = (N/L) × G̃

    where G̃ is the dimensionless 6×6 matrix:

        G̃ = [[0,    0,      0,       0,    0,       0     ],
              [0,    6/5,    L/10,    0,   -6/5,     L/10  ],
              [0,    L/10,   2L²/15,  0,   -L/10,   -L²/30 ],
              [0,    0,      0,       0,    0,       0     ],
              [0,   -6/5,   -L/10,    0,    6/5,    -L/10  ],
              [0,    L/10,  -L²/30,   0,   -L/10,   2L²/15 ]]

    Positive N (tension) increases lateral stiffness.
    Negative N (compression) decreases lateral stiffness → buckling when det(Ke+Kg)=0.

    The local matrix is transformed to global coordinates via T (the standard 6×6
    beam element rotation matrix) and assembled into the ND×ND global DOF space,
    with constrained DOFs omitted.

    Parameters
    ----------
    CORD        (NCT,2) node coordinates
    ECON        (NE,2)  element connectivity (0-indexed)
    DOF         (NCT,3) free DOF flags (1=free, 0=constrained)
    CT          (NE,)   cumulative axial tensions (positive=tension)
    OLEN_elems  (NE,)   original element lengths
    NE, NCT, ND         element / node / DOF counts

    Returns
    -------
    KG  (ND,ND) global geometric stiffness matrix
    """
    KG = np.zeros((ND, ND))

    # Pre-compute the starting global DOF index for each node
    node_dof_start = np.zeros(NCT, dtype=int)
    idx = 0
    for nd in range(NCT):
        node_dof_start[nd] = idx
        idx += int(np.sum(DOF[nd]))

    for el in range(NE):
        N = CT[el]          # axial tension (positive = tension)
        if abs(N) < 1.0e-14:
            continue

        n1, n2 = ECON[el, 0], ECON[el, 1]
        L = OLEN_elems[el]

        # Direction cosines
        dx = CORD[n2, 0] - CORD[n1, 0]
        dy = CORD[n2, 1] - CORD[n1, 1]
        c  = dx / L
        s  = dy / L

        # Local geometric stiffness 6×6
        fac = N / L
        L2  = L * L
        Kg_loc = fac * np.array([
            [ 0,     0,        0,          0,     0,        0       ],
            [ 0,   6/5,     L/10,          0,  -6/5,     L/10      ],
            [ 0,  L/10,  2*L2/15,          0, -L/10,   -L2/30      ],
            [ 0,     0,        0,          0,     0,        0       ],
            [ 0,  -6/5,    -L/10,          0,   6/5,    -L/10      ],
            [ 0,  L/10,   -L2/30,          0, -L/10,   2*L2/15     ]
        ])

        # Coordinate transformation matrix T (6×6)
        T = np.array([
            [ c, -s, 0,  0,  0, 0],
            [ s,  c, 0,  0,  0, 0],
            [ 0,  0, 1,  0,  0, 0],
            [ 0,  0, 0,  c, -s, 0],
            [ 0,  0, 0,  s,  c, 0],
            [ 0,  0, 0,  0,  0, 1]
        ])

        # Global geometric stiffness for this element (6×6)
        Kg_gl = T @ Kg_loc @ T.T

        # Build mapping from local element DOF index (0..5) to global DOF index (-1 if fixed)
        # Local order: [n1_x, n1_y, n1_z,  n2_x, n2_y, n2_z]
        local_to_global = []
        for nd, base in [(n1, node_dof_start[n1]), (n2, node_dof_start[n2])]:
            offset = 0
            for coord in range(3):
                if DOF[nd, coord]:
                    local_to_global.append(base + offset)
                    offset += 1
                else:
                    local_to_global.append(-1)   # constrained — skip

        # Scatter into KG
        for i_loc, gi in enumerate(local_to_global):
            if gi < 0:
                continue
            for j_loc, gj in enumerate(local_to_global):
                if gj < 0:
                    continue
                KG[gi, gj] += Kg_gl[i_loc, j_loc]

    return KG


def solve_with_active_set(KSAT, LV, ND, DOF, RTYPE, NCT,
                          tol_contact=1e-6, max_iter=50, verbose=False):
    """
    Solve with one-way reaction constraints using active set method
    
    Returns:
        disp: displacement vector (ND,)
        active: boolean array indicating active constraints (ND,)
    """
    # Map global DOF indices to node/coord
    dof_to_rtype = []
    
    for nd in range(NCT):
        for coord in range(3):
            if DOF[nd, coord]:
                dof_to_rtype.append(RTYPE[nd, coord])
    
    # Classify DOFs
    free_dofs = [i for i, rt in enumerate(dof_to_rtype) if rt == 0]
    bi_dofs = [i for i, rt in enumerate(dof_to_rtype) if rt == 1]
    pos_dofs = [i for i, rt in enumerate(dof_to_rtype) if rt == 2]
    neg_dofs = [i for i, rt in enumerate(dof_to_rtype) if rt == 3]
    
    # If no one-way reactions, use standard solver
    if len(pos_dofs) == 0 and len(neg_dofs) == 0:
        if verbose:
            print("  No one-way reactions. Using standard solver.")
        active = np.zeros(ND, dtype=bool)
        active[bi_dofs] = True
        
        if len(free_dofs) == 0:
            return np.zeros(ND), active
        
        disp = np.zeros(ND)
        free_dofs = np.array(free_dofs)
        try:
            disp[free_dofs] = solve(KSAT[np.ix_(free_dofs, free_dofs)], LV[free_dofs])
        except (LinAlgError, np.linalg.LinAlgError):
            return None, active
        return disp, active
    
    # Active set iteration
    active_set = set(bi_dofs + pos_dofs + neg_dofs)
    
    if verbose:
        print(f"  Active set iteration:")
        print(f"    Free: {len(free_dofs)}, Bi: {len(bi_dofs)}, "
              f"Pos: {len(pos_dofs)}, Neg: {len(neg_dofs)}")
    
    for iteration in range(max_iter):
        active_list = sorted(list(active_set))
        free_list = sorted(set(range(ND)) - active_set)
        
        if len(free_list) == 0:
            disp = np.zeros(ND)
            Rv = KSAT @ disp - LV
            break
        
        # Solve
        disp = np.zeros(ND)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', LinAlgWarning)
                disp[free_list] = solve(
                    KSAT[np.ix_(free_list, free_list)],
                    LV[free_list] - KSAT[np.ix_(free_list, active_list)] @ disp[active_list]
                )
        except (LinAlgError, np.linalg.LinAlgError):
            active = np.zeros(ND, dtype=bool)
            active[list(active_set)] = True
            return None, active
        
        # Early exit: if displacements are unreasonably large, a mechanism has formed
        if np.max(np.abs(disp)) > 1.0e10:
            if verbose:
                print(f"    Iter {iteration}: collapse mechanism detected (disp={np.max(np.abs(disp)):.2e})")
            active = np.zeros(ND, dtype=bool)
            active[list(active_set)] = True
            return None, active
        
        # Reactions
        Rv = KSAT @ disp - LV
        
        # Check complementarity
        violations = []
        
        for i in pos_dofs:
            if i in active_set:
                if Rv[i] < -tol_contact:
                    violations.append(('remove_pos', i, Rv[i]))
            else:
                if disp[i] < -tol_contact:
                    violations.append(('add_pos', i, disp[i]))
        
        for i in neg_dofs:
            if i in active_set:
                if Rv[i] > tol_contact:
                    violations.append(('remove_neg', i, Rv[i]))
            else:
                if disp[i] > tol_contact:
                    violations.append(('add_neg', i, disp[i]))
        
        max_viol = max([abs(v[2]) for v in violations]) if violations else 0
        
        if verbose and iteration < 10:
            print(f"    Iter {iteration}: {len(active_list)} active, "
                  f"{len(violations)} violations, max={max_viol:.2e}")
        
        if len(violations) == 0:
            if verbose:
                print(f"    Converged in {iteration} iterations")
            break
        
        # Update
        violations.sort(key=lambda x: abs(x[2]), reverse=True)
        vtype, vidx, vval = violations[0]
        
        if 'remove' in vtype:
            active_set.remove(vidx)
        else:
            active_set.add(vidx)
    
    else:
        if verbose:
            print(f"  Warning: Did not converge in {max_iter} iterations")
    
    active = np.zeros(ND, dtype=bool)
    active[list(active_set)] = True
    
    return disp, active

def epframe_oneway_analysis(input_file, output_file, 
                            tol_contact=1e-6, tol_plastic_factor=0.001,
                            max_contact_iter=50):
    """Main elastic-plastic frame analysis with one-way reactions"""
    
    # Read input
    title, NCT, NE, E, Fy, CORD, DOF, RTYPE, ECON, SMA, AREA, ZS, PM, PY, loads = \
        read_input_file(input_file)
    
    # Open output files
    fp = open(output_file, 'w')
    csv_file = output_file + '.csv'
    csv_fp = open(csv_file, 'w')
    
    # Write CSV: title row first, then column headers
    csv_fp.write(f'"{title}"\n')
    csv_fp.write(get_csv_header(NCT, NE) + '\n')
    
    fp.write("%\n")
    fp.write(f"%     {title}\n")
    fp.write("%     " + "-"*len(title) + "\n%\n")
    
    # Calculate number of DOFs
    ND = int(np.sum(DOF))

    # Displacement limit: 10% of the larger frame dimension
    x_range = CORD[:, 0].max() - CORD[:, 0].min()
    y_range = CORD[:, 1].max() - CORD[:, 1].min()
    DLMT = 0.1 * max(x_range, y_range)

    # Degree of static indeterminacy (before any hinges):
    #   DI = (element force DOFs) − (free structural DOFs) = 3·NE − ND
    # Each hinge reduces DI by 1. A mechanism forms when NCYCL ≥ DI + 1.
    # (E3 = 3·NE is computed later; use 3*NE directly here)
    DI = 3 * NE - ND
    
    # Initialize load vector
    LV = np.zeros(ND)
    
    # Process loads
    for load in loads:
        NH, FX, FY, FZ = load
        OLEN = [FX, FY, FZ]
        
        NN = 0
        for nd in range(NH):
            NN += np.sum(DOF[nd])
        
        for k in range(3):
            if DOF[NH, k]:
                LV[NN] = OLEN[k]
                NN += 1
    
    fp.write("%     * GENERAL DATA\n")
    fp.write(f"%          NUMBER OF NODES         {NCT:6d}\n")
    fp.write(f"%          NUMBER OF ELEMENTS      {NE:6d}\n")
    fp.write(f"%          MOD OF ELASTICITY  {E:12.1f}\n")
    fp.write(f"%          YIELD STRESS       {Fy:12.1f}\n")
    fp.write(f"%          STATIC INDETERMINACY    {DI:6d}   "
             f"(mechanism forms after {DI+1} hinges)\n")
    fp.write(f"%          DISPLACEMENT LIMIT {DLMT:12.2f}   "
             f"(0.10 × max frame dimension)\n%\n")
    
    # Input echo
    fp.write("%\n%     * DATA FOR NODES\n")
    fp.write("%           NODE   X-COORD   Y-COORD    RX-TYPE  RY-TYPE  RZ-TYPE\n%\n")
    rtype_names = ['FREE', 'BI', 'POS', 'NEG']
    for i in range(NCT):
        fp.write(f"%{i+1:14d}{CORD[i,0]:12.2f}{CORD[i,1]:10.2f}"
                f"{rtype_names[RTYPE[i,0]]:>10s}{rtype_names[RTYPE[i,1]]:>9s}{rtype_names[RTYPE[i,2]]:>9s}\n")
    
    fp.write("%\n%     * DATA FOR ELEMENTS\n")
    fp.write("%         ELEMENT    N1      N2       IXX      AREA         Z        MP        PY\n%\n")
    for i in range(NE):
        fp.write(f"%{i+1:14d}{ECON[i,0]+1:9d}{ECON[i,1]+1:8d}"
                f"{SMA[i]:10.2f}{AREA[i]:10.2f}{ZS[i]:10.2f}{PM[i]:10.2f}{PY[i]:10.2f}\n")
    
    fp.write("%\n%     * DATA FOR LOADS\n")
    fp.write("%           NODE        PX        PY        PZ\n")
    for load in loads:
        NH, FX, FY, FZ = load
        fp.write(f"%{NH+1:14d}{FX:12.2f}{FY:10.2f}{FZ:10.2f}\n")
    
    # Calculate element lengths
    OLEN_elems = np.zeros(NE)
    for i in range(NE):
        N1, N2 = ECON[i, 0], ECON[i, 1]
        dx = CORD[N1, 0] - CORD[N2, 0]
        dy = CORD[N1, 1] - CORD[N2, 1]
        OLEN_elems[i] = sqrt(dx*dx + dy*dy)
    
    # Calculate stiffness coefficients
    SF = np.zeros((2*NE, 2))
    SA = np.zeros(NE)
    
    for i in range(NE):
        SF[2*i+1, 1] = 4.0 * E * SMA[i] / OLEN_elems[i]
        SF[2*i, 0] = SF[2*i+1, 1]
        SF[2*i+1, 0] = 0.5 * SF[2*i+1, 1]
        SF[2*i, 1] = SF[2*i+1, 0]
        SA[i] = E * AREA[i] / OLEN_elems[i]
    
    E2 = 2 * NE
    E3 = 3 * NE
    
    # Initialize cumulative variables
    CM = np.zeros(E2)
    CT = np.zeros(NE)
    CD = np.zeros(ND)
    
    NCYCL = 0
    CLF = 0.0
    
    # Initial active status
    active = np.ones(ND, dtype=bool)
    for i in range(ND):
        active[i] = False
    dof_idx = 0
    for nd in range(NCT):
        for coord in range(3):
            if DOF[nd, coord]:
                if RTYPE[nd, coord] > 0:
                    active[dof_idx] = True
                dof_idx += 1
    
    # Write initial state
    write_csv_row(csv_fp, 0, 0, 0, 0.0, CD, CM, CT, DOF, RTYPE, active, NCT, NE)
    
    # Build compatibility matrix K
    K = np.zeros((ND, E3))
    
    row_start = 0
    
    for nd in range(NCT):
        for el in range(NE):
            if nd == ECON[el, 0]:
                NF = ECON[el, 1]
                EI_near = 2*el
                EI_far = 2*el + 1
            elif nd == ECON[el, 1]:
                NF = ECON[el, 0]
                EI_near = 2*el + 1
                EI_far = 2*el
            else:
                continue
            
            X = CORD[NF, 0] - CORD[nd, 0]
            Y = CORD[NF, 1] - CORD[nd, 1]
            L = sqrt(X*X + Y*Y)
            S = Y / L
            C = X / L
            axial_col = 2*NE + el
            
            NA = row_start
            if DOF[nd, 0]:
                K[NA, EI_near] = S / L
                K[NA, EI_far] = K[NA, EI_near]
                K[NA, axial_col] = -C
                NA += 1
            
            if DOF[nd, 1]:
                K[NA, EI_near] = -C / L
                K[NA, EI_far] = K[NA, EI_near]
                K[NA, axial_col] = -S
                NA += 1
            
            if DOF[nd, 2]:
                K[NA, EI_near] = 1.0
        
        row_start += int(np.sum(DOF[nd]))
    
    # Main analysis loop
    print("\n=== ELASTIC-PLASTIC ANALYSIS WITH GEOMETRIC STIFFNESS AND ONE-WAY REACTIONS ===\n")
    
    while True:
        NCYCL += 1
        print(f"\n--- Load Increment {NCYCL} ---")
        kg_norm = np.linalg.norm(KG) if 'KG' in dir() else 0.0
        
        # Build element stiffness matrix
        S_full = np.zeros((E3, E3))
        
        for i in range(NE):
            row = 2*i
            S_full[row, row] = SF[2*i, 0]
            S_full[row, row+1] = SF[2*i, 1]
            S_full[row+1, row] = SF[2*i+1, 0]
            S_full[row+1, row+1] = SF[2*i+1, 1]
        
        for i in range(NE):
            S_full[E2+i, E2+i] = SA[i]
        
        # Form elastic stiffness matrix
        KSAT = K @ S_full @ K.T

        # Add geometric stiffness (depends on current cumulative axial forces CT)
        KG = build_geometric_stiffness(CORD, ECON, DOF, CT, OLEN_elems, NE, NCT, ND)
        KSAT += KG

        ke_norm = np.linalg.norm(KSAT - KG)
        kg_norm = np.linalg.norm(KG)
        if ke_norm > 0:
            print(f"  ||Kg|| / ||Ke|| = {kg_norm/ke_norm:.4f}"
                  f"  ({'stiffening' if np.trace(KG) >= 0 else 'softening'}"
                  f", max |CT| = {np.max(np.abs(CT)):.1f} kips)")

        # Check for geometric instability (buckling):
        # A negative minimum eigenvalue of KSAT means the structure has buckled.
        # This is only worth checking when geometric stiffness is non-negligible.
        # We defer this check until after solve_with_active_set returns None so
        # that we can distinguish buckling (negative eigenvalue) from a plastic
        # mechanism (zero eigenvalue due to hinge sequence).
        
        # Solve with active set
        disp, active = solve_with_active_set(
            KSAT, LV, ND, DOF, RTYPE, NCT,
            tol_contact, max_contact_iter, verbose=True
        )
        
        if disp is None:
            # Determine whether failure is geometric instability (buckling) or
            # a plastic mechanism.  Buckling → KSAT has a negative eigenvalue.
            # Mechanism → KSAT is singular (smallest eigenvalue ≈ 0).
            # We check the minimum eigenvalue only when KG is non-trivial.
            is_buckling = False
            if kg_norm > 1.0e-10 * ke_norm:
                eigs = np.linalg.eigvalsh(KSAT)
                if eigs[0] < -1.0e-6 * abs(eigs[-1]):
                    is_buckling = True

            if is_buckling:
                msg = (f"     *** GEOMETRIC INSTABILITY (BUCKLING) IN CYCLE {NCYCL}: "
                       f"min eigenvalue = {eigs[0]:.3e}\n")
            else:
                msg = (f"     *** COLLAPSE MECHANISM DETECTED IN CYCLE {NCYCL}: "
                       f"STIFFNESS MATRIX IS SINGULAR\n")
            fp.write(f"%\n%{msg}%\n")
            print(f"\n{msg.strip()}")
            break
        
        # Calculate forces
        CSAT = K.T @ disp
        SATX = S_full @ CSAT
        
        SATX_e2 = SATX[:E2]   # moment rates (one per element end)
        SATX_ct = SATX[E2:]   # axial force rates (one per element)

        # --- Check for compression yield before finding next hinge ---
        # If any element already has |P| >= Py, stop immediately
        for el in range(NE):
            if abs(CT[el]) >= PY[el]:
                msg = (f"     *** COMPRESSION YIELD IN ELEMENT {el+1}: "
                       f"|P|={abs(CT[el]):.1f} >= Py={PY[el]:.1f}\n")
                fp.write(f"%\n%{msg}%\n")
                print(f"\n{msg.strip()}")
                # Signal clean exit via break-out flag
                SATX = None
                break
        if SATX is None:
            break

        # --- P-M quadratic solve for load factor to next hinge ---
        # At element end k (element el, end 0 or 1):
        #   M0 = CM[k]          cumulative moment
        #   mdot = SATX_e2[k]   moment rate
        #   P0 = CT[el]         cumulative axial force
        #   pdot = SATX_ct[el]  axial force rate
        #
        # Hinge condition: |M0 + α·mdot| = Mp·(1 − ((P0+α·pdot)/Py)²)
        # For the end that is loading toward Mp (TEST = M0·mdot >= 0):
        #   Let s = sign(mdot) if M0≈0, else sign(M0)
        #   s·(M0 + α·mdot) = Mp·(1 − (P0 + α·pdot)²/Py²)
        # Rearranged: A·α² + B·α + C = 0 where
        #   A =  Mp·pdot²/Py²
        #   B =  s·mdot + 2·Mp·P0·pdot/Py²
        #   C =  s·M0 − Mp·(1 − P0²/Py²)
        # When pdot=0 this reduces to the original linear formula α = (Mp−|M0|)/|mdot|

        ALF_pm = np.full(E2, 1.0e10)

        for k in range(E2):
            el = k // 2
            mdot = SATX_e2[k]
            M0   = CM[k]
            pdot = SATX_ct[el]
            P0   = CT[el]
            mp   = PM[el]
            py   = PY[el]

            # Only consider ends where moment is growing toward capacity
            if abs(mdot) < tol_plastic_factor * mp and abs(M0) < mp:
                continue
            if M0 * mdot < 0.0:   # moment relieving, not approaching hinge
                continue

            # Current Mu at this end
            ratio0 = P0 / py
            if abs(ratio0) >= 1.0:
                # Already at or past compression yield — flag handled above
                continue
            Mu0 = mp * (1.0 - ratio0**2)

            if abs(M0) >= Mu0 - tol_plastic_factor * mp:
                # Already at capacity (hinge exists here); skip
                continue

            s = np.sign(M0) if abs(M0) > tol_plastic_factor * mp else np.sign(mdot)

            A = mp * pdot**2 / py**2
            B = s * mdot + 2.0 * mp * P0 * pdot / py**2
            C = s * M0 - mp * (1.0 - (P0/py)**2)

            if abs(A) < 1.0e-14 * abs(B):
                # Linear case (pdot ≈ 0): α = −C / B
                if abs(B) > 1.0e-14:
                    alpha = -C / B
                    if alpha > 0.0:
                        ALF_pm[k] = alpha
            else:
                disc = B**2 - 4.0 * A * C
                if disc < 0.0:
                    continue   # no real root — Mu never reached
                sq = sqrt(disc)
                for alpha in [(-B + sq) / (2.0*A), (-B - sq) / (2.0*A)]:
                    if alpha > 1.0e-10:
                        ALF_pm[k] = min(ALF_pm[k], alpha)

        PHN  = int(np.argmin(ALF_pm))
        SALF = ALF_pm[PHN]

        if SALF > 1.0e9:
            msg = f"     *** NO FURTHER HINGE POSSIBLE IN CYCLE {NCYCL}\n"
            fp.write(f"%\n%{msg}%\n")
            print(f"\n{msg.strip()}")
            break
        
        # Update cumulative values
        SATX *= SALF
        disp *= SALF
        CD += disp
        CLF += SALF
        CM += SATX[:E2]
        CT += SATX[E2:]

        # Check cumulative displacements against frame-size limit
        max_cd = np.max(np.abs(CD))
        if max_cd > DLMT:
            msg = (f"     *** CUMULATIVE DISPLACEMENT {max_cd:.2f} EXCEEDS LIMIT "
                   f"{DLMT:.2f} (0.10 × max frame dimension) IN CYCLE {NCYCL}\n")
            fp.write(f"%\n%{msg}%\n")
            print(f"\n{msg.strip()}")
            break

        # Indeterminacy warning: once hinges exceed DI the structure is
        # kinematically a mechanism — geometric stiffness may mask this.
        if NCYCL > DI:
            msg = (f"     *** WARNING: HINGE COUNT ({NCYCL}) EXCEEDS DEGREE OF "
                   f"STATIC INDETERMINACY ({DI}) — KINEMATIC MECHANISM LIKELY\n")
            fp.write(f"%{msg}")
            print(f"  {msg.strip()}")
        
        # Hinge location
        EL = PHN // 2
        if PHN % 2 == 0:
            NH = ECON[EL, 0]
        else:
            NH = ECON[EL, 1]
        
        EL_hinge = EL
        NH_hinge = NH
        
        # Output
        fp.write("%\n%\n%\n")
        fp.write(f"%     * PLASTIC HINGE {NCYCL:3d} FORMED IN ELEMENT {EL+1:3d} "
                f"NEAR NODE {NH+1:3d} WHEN LOAD FACTOR IS {CLF:12.3f}\n%\n")
        
        print(f"*** HINGE {NCYCL} in element {EL+1} near node {NH+1} at λ = {CLF:.3f}")
        
        # Active support status
        fp.write("%     ACTIVE SUPPORT STATUS:\n")
        print("\nActive Support Status:")
        
        dof_idx = 0
        rtype_names = ['FREE', 'BI', 'POS', 'NEG']
        coord_names = ['X', 'Y', 'Z']
        
        for nd in range(NCT):
            has_reaction = np.any(RTYPE[nd, :] > 0)
            if not has_reaction:
                for coord in range(3):
                    if DOF[nd, coord]:
                        dof_idx += 1
                continue
            
            status_parts = []
            for coord in range(3):
                if DOF[nd, coord]:
                    if RTYPE[nd, coord] > 0:
                        is_active = active[dof_idx]
                        status = "ACTIVE" if is_active else "LIFT-OFF"
                        status_parts.append(f"{coord_names[coord]}={rtype_names[RTYPE[nd,coord]]}:{status}")
                    dof_idx += 1
            
            if status_parts:
                status_str = ', '.join(status_parts)
                fp.write(f"%       NODE {nd+1:3d}: {status_str}\n")
                print(f"  Node {nd+1}: {status_str}")
        
        fp.write("%\n")
        
        # Displacements
        fp.write("%          CUMULATIVE DEFORMATIONS\n")
        fp.write("%                NODE    X-DISP       Y-DISP       ROTN\n")
        
        NN = 0
        for nd in range(NCT):
            disp_out = [0.0, 0.0, 0.0]
            for dof in range(3):
                if DOF[nd, dof]:
                    # Negate rotation (dof=2): convert internal CW-positive to CCW-positive
                    disp_out[dof] = -CD[NN] if dof == 2 else CD[NN]
                    NN += 1
            fp.write(f"%{nd+1:19d}{disp_out[0]:13.5f}{disp_out[1]:13.5f}{disp_out[2]:13.5f}\n")
        
        # Moments — include Mu and P/Py
        fp.write("%\n%          CUMULATIVE MOMENTS\n")
        fp.write("%             ELEMENT       END MOMENTS             NODES"
                 "     PLASTIC MOM   ULTIMATE MOM      P / Py\n")
        
        for el in range(NE):
            k = 2*el
            p_ratio = CT[el] / PY[el]
            Mu = PM[el] * (1.0 - p_ratio**2)
            fp.write(f"%{el+1:19d}{-CM[k]:14.2f}{-CM[k+1]:11.2f}"
                    f"{ECON[el,0]+1:7d} AND {ECON[el,1]+1:2d}"
                    f"{PM[el]:14.2f}{Mu:14.2f}{p_ratio:11.4f}\n")
        
        # Axial forces
        fp.write("%\n%          CUMULATIVE TENSION FORCES\n")
        fp.write("%             ELEMENT     TENSION\n")
        
        for el in range(NE):
            fp.write(f"%{el+1:19d}{CT[el]:15.2f}\n")
        
        # Reactions
        fp.write("%\n%          REACTIONS AT SUPPORTS\n")
        fp.write("%                NODE       FX           FY           MZ         STATUS\n")
        
        dof_idx = 0
        for nd in range(NCT):
            is_support = np.any(RTYPE[nd, :] > 0)
            if not is_support:
                for coord in range(3):
                    if DOF[nd, coord]:
                        dof_idx += 1
                continue
            
            Rx, Ry, Rz = 0.0, 0.0, 0.0
            
            for el in range(NE):
                N1, N2 = ECON[el, 0], ECON[el, 1]
                
                if N1 != nd and N2 != nd:
                    continue
                
                dx = CORD[N2, 0] - CORD[N1, 0]
                dy = CORD[N2, 1] - CORD[N1, 1]
                L_el = sqrt(dx*dx + dy*dy)
                Cx = dx / L_el
                Sy = dy / L_el
                
                Me1 = CM[2*el]
                Me2 = CM[2*el + 1]
                N_ax = CT[el]
                V = (Me1 + Me2) / L_el
                
                if N1 == nd:
                    Rx += N_ax * Cx - V * Sy
                    Ry += N_ax * Sy + V * Cx
                    Rz += Me1
                else:
                    Rx += -N_ax * Cx + V * Sy
                    Ry += -N_ax * Sy - V * Cx
                    Rz += Me2
            
            # Check inactive constraints
            inactive = []
            for coord in range(3):
                if DOF[nd, coord] and RTYPE[nd, coord] > 0:
                    if not active[dof_idx]:
                        inactive.append(coord_names[coord])
                    dof_idx += 1
                elif not DOF[nd, coord]:
                    pass
                else:
                    dof_idx += 1
            
            status_str = "LIFT-OFF:" + ','.join(inactive) if inactive else "ACTIVE"
            fp.write(f"%{nd+1:19d}{-Rx:13.2f}{-Ry:13.2f}{Rz:13.2f}    {status_str}\n")
        
        # CSV output
        write_csv_row(csv_fp, NCYCL, EL_hinge+1, NH_hinge+1, CLF, CD, CM, CT, 
                     DOF, RTYPE, active, NCT, NE)
        
        # Modify stiffness
        if PHN % 2 == 0:
            SF[PHN+1, 1] = 0.75 * SF[PHN+1, 1]
            SF[PHN+1, 0] = 0.0
            SF[PHN, 0] = 0.0
            SF[PHN, 1] = 0.0
        else:
            SF[PHN-1, 0] = 0.75 * SF[PHN-1, 0]
            SF[PHN-1, 1] = 0.0
            SF[PHN, 0] = 0.0
            SF[PHN, 1] = 0.0
        
        print(f"    LOAD FACTOR {NCYCL:2d} = {CLF:.3f}")
    
    fp.write(f"%\n%     ANALYSIS COMPLETED: {title} AT LOAD FACTOR {CLF:.3f}\n\n")
    fp.close()
    csv_fp.close()
    
    print(f"\nAnalysis completed at load factor {CLF:.3f}")
    print(f"Results written to {output_file}")
    print(f"Compact data written to {csv_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python epframe_oneway.py input_file output_file")
        sys.exit(1)
    
    epframe_oneway_analysis(sys.argv[1], sys.argv[2])
