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
from math import sqrt
from scipy.linalg import solve

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
    """Write a data row to CSV file"""
    values = [f'{NCYCL}', f'{EL}', f'{NH}', f'{CLF:.6E}']
    
    # CD - cumulative displacement
    NN = 0
    for nd in range(NCT):
        for dof in range(3):
            if DOF[nd, dof]:
                values.append(f'{CD[NN]:.6E}')
                NN += 1
            else:
                values.append(f'{0.0:.6E}')
    
    # CM values
    for i in range(2*NE):
        values.append(f'{CM[i]:.6E}')
    
    # CT values
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
        'X Y Z' or 'XYZ' or 'X* Y* Z*' -> DOF=[0,0,0], RTYPE=[1,1,1]
        'X Y *' or 'XY*' -> DOF=[0,0,1], RTYPE=[1,1,0]
        'X+ Y Z' -> DOF=[0,0,0], RTYPE=[2,1,1]
        '0 0 0' or '* * *' -> DOF=[1,1,1], RTYPE=[0,0,0]
    
    Returns:
        DOF: [DFX, DFY, DFZ] where 1=free, 0=has some constraint
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
        
        # Determine which coordinate this applies to
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
            # Assume numeric old format: 0 or 1
            try:
                val = int(part)
                # This is a position-based value (old format compatibility)
                continue
            except ValueError:
                continue
        
        if coord_idx is None:
            continue
        
        # Parse the constraint character
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
    
    # Invert DOF logic: if we have a reaction type > 0, DOF should be 0
    # DOF indicates if it's a free coordinate (1) or constrained (0)
    for i in range(3):
        if RTYPE[i] > 0:
            DOF[i] = 0
        else:
            DOF[i] = 1
    
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
    FN = int(lines[idx].split()[0])
    idx += 1
    
    parts = lines[idx].split()
    NCT, NE, E = int(parts[0]), int(parts[1]), float(parts[2])
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
        
        # Check if old numeric format or new string format
        if len(parts) >= 6 and parts[3].replace('-','').replace('.','').isdigit():
            # Old numeric format: convert
            DOF[node_num, 0] = 1 - int(parts[3])  # Invert: RCT 1 -> DOF 0
            DOF[node_num, 1] = 1 - int(parts[4])
            DOF[node_num, 2] = 1 - int(parts[5])
            # Assume bidirectional for old format
            RTYPE[node_num, :] = [1 if d == 0 else 0 for d in DOF[node_num, :]]
        else:
            # New format: parse fixity string
            fixity_str = ' '.join(parts[3:])
            dof, rtype = parse_fixity(fixity_str)
            DOF[node_num, :] = dof
            RTYPE[node_num, :] = rtype
        
        idx += 1
    
    # Element data
    ECON = np.zeros((NE, 2), dtype=int)
    SMA = np.zeros(NE)
    AREA = np.zeros(NE)
    PM = np.zeros(NE)
    
    for i in range(NE):
        parts = lines[idx].split()
        elem_num = int(parts[0]) - 1
        ECON[elem_num, 0] = int(parts[1]) - 1
        ECON[elem_num, 1] = int(parts[2]) - 1
        SMA[elem_num] = float(parts[3])
        AREA[elem_num] = float(parts[4])
        PM[elem_num] = float(parts[5])
        idx += 1
    
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
    
    return FN, NCT, NE, E, CORD, DOF, RTYPE, ECON, SMA, AREA, PM, loads

def solve_with_active_set(KSAT, LV, ND, DOF, RTYPE, NCT,
                          tol_contact=1e-6, max_iter=50, verbose=False):
    """
    Solve with one-way reaction constraints using active set method
    
    Returns:
        disp: displacement vector (ND,)
        active: boolean array indicating active constraints (ND,)
    """
    # Map global DOF indices to node/coord
    dof_to_node = []
    dof_to_coord = []
    dof_to_rtype = []
    
    for nd in range(NCT):
        for coord in range(3):
            if DOF[nd, coord]:
                dof_to_node.append(nd)
                dof_to_coord.append(coord)
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
        disp[free_dofs] = solve(KSAT[np.ix_(free_dofs, free_dofs)], LV[free_dofs])
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
        disp[free_list] = solve(
            KSAT[np.ix_(free_list, free_list)],
            LV[free_list] - KSAT[np.ix_(free_list, active_list)] @ disp[active_list]
        )
        
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
    FN, NCT, NE, E, CORD, DOF, RTYPE, ECON, SMA, AREA, PM, loads = \
        read_input_file(input_file)
    
    # Open output files
    fp = open(output_file, 'w')
    csv_file = output_file + '.csv'
    csv_fp = open(csv_file, 'w')
    
    # Write CSV header
    csv_fp.write(get_csv_header(NCT, NE) + '\n')
    
    fp.write("%\n")
    fp.write(f"%     ELASTIC PLASTIC ANALYSIS WITH ONE-WAY REACTIONS - FRAME NO {FN}\n")
    fp.write("%     ---------------------------------------------------------------\n%\n")
    
    # Calculate number of DOFs
    ND = int(np.sum(DOF))
    
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
    fp.write(f"%          MOD OF ELASTICITY  {E:12.1f}\n%\n")
    
    # Input echo with reaction types
    fp.write("%\n%     * DATA FOR NODES\n")
    fp.write("%           NODE   X-COORD   Y-COORD    RX-TYPE  RY-TYPE  RZ-TYPE\n%\n")
    rtype_names = ['FREE', 'BI', 'POS', 'NEG']
    for i in range(NCT):
        fp.write(f"%{i+1:14d}{CORD[i,0]:12.2f}{CORD[i,1]:10.2f}"
                f"{rtype_names[RTYPE[i,0]]:>10s}{rtype_names[RTYPE[i,1]]:>9s}{rtype_names[RTYPE[i,2]]:>9s}\n")
    
    # Element echo
    fp.write("%\n%     * DATA FOR ELEMENTS\n")
    fp.write("%         ELEMENT    N1      N2       IXX      AREA        MP\n%\n")
    for i in range(NE):
        fp.write(f"%{i+1:14d}{ECON[i,0]+1:9d}{ECON[i,1]+1:8d}"
                f"{SMA[i]:10.2f}{AREA[i]:10.2f}{PM[i]:10.2f}\n")
    
    # Load echo
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
    
    # Initial active status (all reactions active)
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
    print("\n=== ELASTIC-PLASTIC ANALYSIS WITH ONE-WAY REACTIONS ===\n")
    
    while True:
        NCYCL += 1
        print(f"\n--- Load Increment {NCYCL} ---")
        
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
        
        # Form stiffness matrix
        KSAT = K @ S_full @ K.T
        
        # Solve with active set for one-way reactions
        disp, active = solve_with_active_set(
            KSAT, LV, ND, DOF, RTYPE, NCT,
            tol_contact, max_contact_iter, verbose=True
        )
        
        # Check deformations
        DLMT = 1000.0
        max_disp = np.max(np.abs(disp))
        
        if max_disp > DLMT:
            fp.write(f"%\n%     *** DEFORMATIONS LARGER THAN {DLMT:.1f} IN CYCLE NO {NCYCL:4d}\n%\n")
            print(f"\nDeformations exceed limit ({max_disp:.2e} > {DLMT:.1f})")
            break
        
        # Calculate forces
        CSAT = K.T @ disp
        SATX = S_full @ CSAT
        
        # Load factors to plastic hinge
        SATX_e2 = SATX[:E2]
        elem_indices = np.arange(E2) // 2
        PM_expanded = PM[elem_indices]
        
        ZERO = tol_plastic_factor * PM_expanded
        small_satx = np.abs(SATX_e2) < ZERO
        safe_denom = np.where(small_satx, 1.0, np.abs(SATX_e2))
        ALF = np.where(small_satx, 1.0E10, (PM_expanded - np.abs(CM)) / safe_denom)
        
        TEST = CM * SATX_e2
        valid_mask = TEST >= 0.0
        ALF_masked = np.where(valid_mask, ALF, 1.0E10)
        PHN = int(np.argmin(ALF_masked))
        SALF = ALF_masked[PHN]
        
        # Update cumulative values
        SATX *= SALF
        disp *= SALF
        CD += disp
        CLF += SALF
        CM += SATX[:E2]
        CT += SATX[E2:]
        
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
                    disp_out[dof] = CD[NN]
                    NN += 1
            fp.write(f"%{nd+1:19d}{disp_out[0]:13.5f}{disp_out[1]:13.5f}{disp_out[2]:13.5f}\n")
        
        # Moments
        fp.write("%\n%          CUMULATIVE MOMENTS\n")
        fp.write("%             ELEMENT       END MOMENTS             NODES     PLASTIC MOM\n")
        
        for el in range(NE):
            k = 2*el
            fp.write(f"%{el+1:19d}{CM[k]:14.2f}{CM[k+1]:11.2f}"
                    f"{ECON[el,0]+1:7d} AND{ECON[el,1]+1:2d}{PM[el]:14.2f}\n")
        
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
            
            # Check if any constraints are inactive
            inactive = []
            for coord in range(3):
                if DOF[nd, coord] and RTYPE[nd, coord] > 0:
                    if not active[dof_idx]:
                        inactive.append(coord_names[coord])
                    dof_idx += 1
                elif not DOF[nd, coord]:
                    pass  # Fixed DOF, not tracked
                else:
                    dof_idx += 1
            
            status_str = "LIFT-OFF:" + ','.join(inactive) if inactive else "ACTIVE"
            fp.write(f"%{nd+1:19d}{-Rx:13.2f}{-Ry:13.2f}{-Rz:13.2f}    {status_str}\n")
        
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
    
    fp.write(f"%\n%     ANALYSIS COMPLETED FOR FRAME NO {FN:3d} AT LOAD FACTOR {CLF:.3f}\n\n")
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
"""
ELASTIC PLASTIC ANALYSIS OF A PLANE FRAME WITH ONE-WAY REACTIONS
Combines plastic hinge formation with contact constraint handling via QP

One-way reaction types:
  '0' or 0   : No reaction (free)
  '*' or 1   : Bi-directional reaction (standard support)
  '+' or 2   : Reaction only in positive direction
  '-' or 3   : Reaction only in negative direction
"""

import numpy as np
import sys
from math import sqrt
from scipy.linalg import solve

class Index1Array:
    """Wrapper to provide 1-based indexing for numpy arrays"""
    def __init__(self, shape, dtype=float):
        if isinstance(shape, (int, np.integer)):
            self._data = np.zeros(shape + 1, dtype=dtype)
            self._offset = 1
        elif isinstance(shape, (tuple, list)):
            if len(shape) == 1:
                self._data = np.zeros(shape[0] + 1, dtype=dtype)
                self._offset = 1
            elif len(shape) == 2:
                self._data = np.zeros((shape[0] + 1, shape[1] + 1), dtype=dtype)
                self._offset = (1, 1)
        else:
            raise ValueError("Shape must be int or tuple/list")
    
    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self._data[key[0], key[1]]
        return self._data[key]
    
    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            self._data[key[0], key[1]] = value
        else:
            self._data[key] = value

def parse_fixity(fixity_str):
    """
    Parse fixity string to DOF flags and reaction type
    Returns: (DFX, DFY, DFZ, RX_type, RY_type, RZ_type)
      DF: 1=free, 0=fixed (for compatibility matrix)
      R_type: 0=free, 1=bidirectional, 2=positive, 3=negative
    """
    fixity_str = str(fixity_str).upper().strip()
    
    # Map character to reaction type
    def char_to_type(c):
        if c == '0':
            return 1, 0  # free: DF=1, type=0
        elif c == '*' or c == '1':
            return 0, 1  # bidirectional: DF=0, type=1
        elif c == '+' or c == '2':
            return 0, 2  # positive only: DF=0, type=2
        elif c == '-' or c == '3':
            return 0, 3  # negative only: DF=0, type=3
        else:
            raise ValueError(f"Invalid fixity character: {c}")
    
    parts = fixity_str.split()
    
    # Initialize all as free
    DF = [1, 1, 1]
    R_type = [0, 0, 0]
    
    # Process each part
    for part in parts:
        if 'X' in part:
            DF[0], R_type[0] = char_to_type(part.replace('X', ''))
        elif 'Y' in part:
            DF[1], R_type[1] = char_to_type(part.replace('Y', ''))
        elif 'Z' in part:
            DF[2], R_type[2] = char_to_type(part.replace('Z', ''))
    
    return DF[0], DF[1], DF[2], R_type[0], R_type[1], R_type[2]

def read_input_file(filename):
    """Read input file with comment support and one-way reactions"""
    with open(filename, 'r') as f:
        lines = []
        for line in f:
            line = line.split('#')[0].strip()
            if line:
                lines.append(line)
    
    idx = 0
    JFN = int(lines[idx].split()[0])
    idx += 1
    
    parts = lines[idx].split()
    JCT, NM, E = int(parts[0]), int(parts[1]), float(parts[2])
    idx += 1
    
    # Node data
    CORD = Index1Array((JCT, 2))
    JTYPE = Index1Array((JCT, 3), dtype=int)
    RTYPE = Index1Array((JCT, 3), dtype=int)  # Reaction types
    
    for i in range(1, JCT+1):
        parts = lines[idx].split()
        node_num = int(parts[0])
        CORD[i,1] = float(parts[1])
        CORD[i,2] = float(parts[2])
        
        # Parse fixity string
        fixity_str = ' '.join(parts[3:])
        dfx, dfy, dfz, rtx, rty, rtz = parse_fixity(fixity_str)
        
        JTYPE[i,1] = dfx
        JTYPE[i,2] = dfy
        JTYPE[i,3] = dfz
        
        RTYPE[i,1] = rtx
        RTYPE[i,2] = rty
        RTYPE[i,3] = rtz
        
        idx += 1
    
    # Member data
    MCON = Index1Array((NM, 2), dtype=int)
    SMA = Index1Array(NM)
    AREA = Index1Array(NM)
    PM = Index1Array(NM)
    
    for i in range(1, NM+1):
        parts = lines[idx].split()
        mem_num = int(parts[0])
        MCON[i,1] = int(parts[1])
        MCON[i,2] = int(parts[2])
        SMA[i] = float(parts[3])
        AREA[i] = float(parts[4])
        PM[i] = float(parts[5])
        idx += 1
    
    # Load data
    LN = int(lines[idx].split()[0])
    idx += 1
    
    loads = []
    for i in range(LN):
        parts = lines[idx].split()
        JN = int(parts[0])
        FX = float(parts[1])
        FY = float(parts[2])
        FZ = float(parts[3])
        loads.append((JN, FX, FY, FZ))
        idx += 1
    
    return JFN, JCT, NM, E, CORD, JTYPE, RTYPE, MCON, SMA, AREA, PM, loads

def solve_with_active_set(K_np, VL_np, L, JTYPE, RTYPE, JCT, 
                          tol_contact=1e-6, max_iter=50, verbose=True):
    """
    Solve displacement problem with one-way reaction constraints
    
    Returns:
    --------
    disp_np : ndarray (L,)
        Displacements at free DOFs
    active : ndarray (L,)
        Boolean array of active constraints
    """
    # Identify constraint types for each DOF
    free_dofs = []
    bi_dofs = []
    pos_dofs = []
    neg_dofs = []
    
    dof_idx = 0
    for i in range(1, JCT+1):
        for j in range(1, 4):
            if JTYPE[i,j]:  # This is a free DOF
                rtype = RTYPE[i,j]
                if rtype == 0:
                    free_dofs.append(dof_idx)
                elif rtype == 1:
                    bi_dofs.append(dof_idx)
                elif rtype == 2:
                    pos_dofs.append(dof_idx)
                elif rtype == 3:
                    neg_dofs.append(dof_idx)
                dof_idx += 1
    
    free_dofs = np.array(free_dofs)
    bi_dofs = np.array(bi_dofs)
    pos_dofs = np.array(pos_dofs)
    neg_dofs = np.array(neg_dofs)
    
    # If no one-way reactions, use standard solver
    if len(pos_dofs) == 0 and len(neg_dofs) == 0:
        if verbose:
            print("  No one-way reactions. Using standard solver.")
        
        # All constraints are active
        active = np.zeros(L, dtype=bool)
        active[bi_dofs] = True
        
        if len(free_dofs) == 0:
            return np.zeros(L), active
        
        disp_np = np.zeros(L)
        disp_np[free_dofs] = solve(K_np[np.ix_(free_dofs, free_dofs)], 
                                    VL_np[free_dofs])
        return disp_np, active
    
    # Active set iteration for one-way reactions
    active_set = set(bi_dofs) | set(pos_dofs) | set(neg_dofs)
    
    if verbose:
        print(f"  Starting active set iteration:")
        print(f"    Free DOFs: {len(free_dofs)}")
        print(f"    Bidirectional: {len(bi_dofs)}")
        print(f"    Positive-only: {len(pos_dofs)}")
        print(f"    Negative-only: {len(neg_dofs)}")
    
    for iteration in range(max_iter):
        active_list = sorted(list(active_set))
        free_list = sorted(set(range(L)) - active_set)
        
        if len(free_list) == 0:
            disp_np = np.zeros(L)
            Rv_np = K_np @ disp_np - VL_np
            break
        
        # Solve for displacements
        disp_np = np.zeros(L)
        disp_np[free_list] = solve(
            K_np[np.ix_(free_list, free_list)],
            VL_np[free_list] - K_np[np.ix_(free_list, active_list)] @ disp_np[active_list]
        )
        
        # Calculate reactions
        Rv_np = K_np @ disp_np - VL_np
        
        # Check complementarity
        violations = []
        
        for i in pos_dofs:
            if i in active_set:
                if Rv_np[i] < -tol_contact:
                    violations.append(('remove_pos', i, Rv_np[i]))
            else:
                if disp_np[i] < -tol_contact:
                    violations.append(('add_pos', i, disp_np[i]))
        
        for i in neg_dofs:
            if i in active_set:
                if Rv_np[i] > tol_contact:
                    violations.append(('remove_neg', i, Rv_np[i]))
            else:
                if disp_np[i] > tol_contact:
                    violations.append(('add_neg', i, disp_np[i]))
        
        max_viol = max([abs(v[2]) for v in violations]) if violations else 0
        
        if verbose:
            print(f"    Iter {iteration}: {len(active_list)} active, "
                  f"{len(violations)} violations, max={max_viol:.2e}")
        
        if len(violations) == 0:
            if verbose:
                print(f"    Converged in {iteration} iterations")
            break
        
        # Update active set
        violations.sort(key=lambda x: abs(x[2]), reverse=True)
        vtype, vidx, vval = violations[0]
        
        if 'remove' in vtype:
            active_set.remove(vidx)
        else:
            active_set.add(vidx)
    
    else:
        print(f"  Warning: Active set did not converge in {max_iter} iterations")
    
    # Create active array
    active = np.zeros(L, dtype=bool)
    active[list(active_set)] = True
    
    return disp_np, active

def epframe_oneway_analysis(input_file, output_file, tol_contact=1e-6, 
                            tol_plastic_factor=0.001, max_contact_iter=50):
    """
    Main elastic-plastic frame analysis with one-way reactions
    """
    # Read input
    JFN, JCT, NM, E, CORD, JTYPE, RTYPE, MCON, SMA, AREA, PM, loads = \
        read_input_file(input_file)
    
    # Open output file
    fp = open(output_file, 'w')
    
    fp.write("%\n")
    fp.write(f"%     ELASTIC PLASTIC ANALYSIS WITH ONE-WAY REACTIONS - FRAME NO {JFN:3d}\n")
    fp.write("%     ----------------------------------------------------------------\n%\n")
    
    # Calculate number of DOFs
    L = 0
    for i in range(1, JCT+1):
        for j in range(1, 4):
            L += JTYPE[i,j]
    
    # Initialize load vector
    VL = Index1Array(L)
    for i in range(1, L+1):
        VL[i] = 0.0
    
    # Process loads
    OLEN = Index1Array(3)
    for load in loads:
        JN, FX, FY, FZ = load
        OLEN[1], OLEN[2], OLEN[3] = FX, FY, FZ
        
        LL = 0
        LJ = JN - 1
        if LJ > 0:
            for j in range(1, LJ+1):
                for k in range(1, 4):
                    LL += JTYPE[j,k]
        
        for k in range(1, 4):
            if JTYPE[JN,k]:
                LL += 1
                VL[LL] = OLEN[k]
    
    fp.write("%     * GENERAL DATA\n")
    fp.write(f"%          NUMBER OF NODES           {JCT:6d}\n")
    fp.write(f"%          NUMBER OF MEMBERS         {NM:6d}\n")
    fp.write(f"%          MOD OF ELASTICITY    {E:12.1f}\n%\n")
    
    # Calculate member lengths
    OLEN_members = Index1Array(NM)
    for i in range(1, NM+1):
        J1 = MCON[i,1]
        J2 = MCON[i,2]
        X = CORD[J1,1] - CORD[J2,1]
        Y = CORD[J1,2] - CORD[J2,2]
        OLEN_members[i] = sqrt(X*X + Y*Y)
    
    # Calculate stiffness coefficients
    SF = Index1Array((2*NM, 2))
    SA = Index1Array(NM)
    
    for i in range(1, NM+1):
        SF[2*i, 2] = 4.0 * E * SMA[i] / OLEN_members[i]
        SF[2*i-1, 1] = SF[2*i, 2]
        SF[2*i, 1] = 0.5 * SF[2*i, 2]
        SF[2*i-1, 2] = SF[2*i, 1]
        SA[i] = E * AREA[i] / OLEN_members[i]
    
    M2 = 2 * NM
    M3 = 3 * NM
    
    # Initialize cumulative variables
    CM = Index1Array(M2)
    CT = Index1Array(NM)
    CX = Index1Array(L)
    for i in range(1, M2+1):
        CM[i] = 0.0
    for i in range(1, NM+1):
        CT[i] = 0.0
    for i in range(1, L+1):
        CX[i] = 0.0
    
    NCYCL = 0
    CLG = 0.0
    
    # Build compatibility matrix K
    K = Index1Array((L, M3))
    for i in range(1, L+1):
        for j in range(1, M3+1):
            K[i,j] = 0.0
    
    NJ = 0
    NK = 0
    
    for J in range(1, JCT+1):
        for M in range(1, NM+1):
            NA = NJ
            
            if J == MCON[M,1]:
                JF = MCON[M,2]
                MJ = 2*M - 1
                MF = MJ + 1
            elif J == MCON[M,2]:
                JF = MCON[M,1]
                MJ = 2*M
                MF = MJ - 1
            else:
                continue
            
            X = CORD[JF,1] - CORD[J,1]
            Y = CORD[JF,2] - CORD[J,2]
            D = sqrt(X*X + Y*Y)
            S = Y / D
            C = X / D
            NN = 2*NM + M
            
            if JTYPE[J,1]:
                NA = NA + 1
                K[NA,MJ] = S / D
                K[NA,MF] = K[NA,MJ]
                K[NA,NN] = -C
            
            if JTYPE[J,2]:
                NA = NA + 1
                K[NA,MJ] = -C / D
                K[NA,MF] = K[NA,MJ]
                K[NA,NN] = -S
            
            if JTYPE[J,3]:
                NA = NA + 1
                K[NA,MJ] = 1.0
            
            if NA > NK:
                NK = NA
        
        NJ = NK
    
    # Convert to numpy for efficient computation
    K_np = K._data[1:L+1, 1:M3+1]
    VL_np = VL._data[1:L+1]
    
    # Main analysis loop
    print("\n=== STARTING ELASTIC-PLASTIC ANALYSIS WITH ONE-WAY REACTIONS ===\n")
    
    while True:
        NCYCL += 1
        print(f"\n--- Load Increment {NCYCL} ---")
        
        # Form stiffness matrix using vectorized operations
        SK = np.zeros((L, M3))
        SF_np = SF._data[1:M2+1, 1:3]
        SA_np = SA._data[1:NM+1]
        
        for j in range(L):
            for i in range(M2):
                k = ((i+1)//2)*2 - 1
                SK[j, i] = SF_np[i, 0]*K_np[j, k-1] + SF_np[i, 1]*K_np[j, k]
            for i in range(NM):
                SK[j, M2+i] = SA_np[i] * K_np[j, M2+i]
        
        # Form KSAT = K^T @ SK
        KSAT_np = K_np.T @ SK
        
        # Solve with active set for one-way reactions
        disp_np, active = solve_with_active_set(
            KSAT_np, VL_np, L, JTYPE, RTYPE, JCT,
            tol_contact, max_contact_iter, verbose=True
        )
        
        # Check for excessive deformations
        XLMT = 1000.0
        max_disp = np.max(np.abs(disp_np))
        
        if max_disp > XLMT:
            fp.write(f"%\n%     *** DEFORMATIONS LARGER THAN {XLMT:.1f} IN CYCLE NO {NCYCL:4d}\n%\n")
            print(f"\nDeformations exceed limit ({max_disp:.2e} > {XLMT:.1f})")
            break
        
        # Calculate member forces
        CSAT_np = K_np.T @ disp_np
        
        # Calculate moments and axial forces
        SATX = Index1Array(M3)
        for I in range(1, M2+1):
            k = ((I+1)//2)*2 - 1
            SATX[I] = SF[I,1]*CSAT_np[k-1] + SF[I,2]*CSAT_np[k]
        
        for I in range(1, NM+1):
            k = M2 + I
            SATX[k] = SA[I] * CSAT_np[k-1]
        
        # Calculate load factors to plastic hinge
        ALG = Index1Array(M2)
        PM_np = PM._data[1:NM+1]
        CM_np = CM._data[1:M2+1]
        SATX_np = SATX._data[1:M2+1]
        
        for I in range(1, M2+1):
            k = (I+1)//2
            ZERO = tol_plastic_factor * PM[k]
            if abs(SATX[I]) < ZERO:
                ALG[I] = 1.0E10
            else:
                ALG[I] = (PM[k] - abs(CM[I])) / abs(SATX[I])
        
        # Find minimum load factor
        SALG = 1.0E10
        NPH = 0
        TEST_np = CM_np * SATX_np
        ALG_np = ALG._data[1:M2+1]
        
        valid_mask = TEST_np >= 0.0
        if np.any(valid_mask):
            valid_ALG = np.where(valid_mask, ALG_np, 1.0E10)
            min_idx = np.argmin(valid_ALG)
            SALG = valid_ALG[min_idx]
            NPH = min_idx + 1
        
        # Scale forces and displacements
        SATX_np = SATX._data[1:M3+1] * SALG
        for I in range(1, M3+1):
            SATX[I] = SATX_np[I-1]
        
        CLG += SALG
        
        # Update cumulative values
        CM_np = CM._data[1:M2+1] + SATX_np[:M2]
        for I in range(1, M2+1):
            CM[I] = CM_np[I-1]
        
        CT_np = CT._data[1:NM+1] + SATX_np[M2:M2+NM]
        for I in range(1, NM+1):
            CT[I] = CT_np[I-1]
        
        disp_scaled_np = disp_np * SALG
        CX_np = CX._data[1:L+1] + disp_scaled_np
        for I in range(1, L+1):
            CX[I] = CX_np[I-1]
        
        # Determine hinge location
        I = (NPH + 1) // 2
        k = (NPH // 2) * 2 - NPH
        if k:
            J = MCON[I,1]
        else:
            J = MCON[I,2]
        
        II = I
        JJ = J
        
        # Print results
        fp.write("%\n%\n%\n")
        fp.write(f"%     * PLASTIC HINGE {NCYCL:3d} FORMED IN MEMBER {I:3d} NEAR NODE {J:3d} WHEN LOAD FACTOR IS {CLG:12.3f}\n%\n")
        
        print(f"\n*** PLASTIC HINGE {NCYCL} formed in member {I} near node {J} at λ = {CLG:.3f}")
        
        # Print active support status
        fp.write("%     ACTIVE SUPPORT STATUS:\n")
        print("\nActive Support Status:")
        
        dof_idx = 0
        for i in range(1, JCT+1):
            has_reaction = False
            status_str = []
            
            for j in range(1, 4):
                if JTYPE[i,j] == 0:  # Fixed DOF
                    continue
                
                rtype = RTYPE[i,j]
                if rtype > 0:  # Some type of reaction
                    has_reaction = True
                    coord = ['X', 'Y', 'Z'][j-1]
                    is_active = active[dof_idx]
                    
                    if rtype == 1:
                        type_str = "BIDIRECTIONAL"
                    elif rtype == 2:
                        type_str = "POSITIVE-ONLY"
                    elif rtype == 3:
                        type_str = "NEGATIVE-ONLY"
                    
                    status = "ACTIVE" if is_active else "LIFT-OFF"
                    status_str.append(f"{coord}={type_str}:{status}")
                
                dof_idx += 1
            
            if has_reaction:
                fp.write(f"%       NODE {i:3d}: {', '.join(status_str)}\n")
                print(f"  Node {i}: {', '.join(status_str)}")
        
        fp.write("%\n")
        
        # Print displacements
        fp.write("%          CUMULATIVE DEFORMATIONS\n")
        fp.write("%               NODE    X-DISP       Y-DISP       ROTN\n")
        
        LL = 0
        for I in range(1, JCT+1):
            disp_out = [0.0, 0.0, 0.0]
            for J in range(1, 4):
                if JTYPE[I,J]:
                    LL += 1
                    disp_out[J-1] = CX[LL]
            fp.write(f"%{I:19d}{disp_out[0]:13.5f}{disp_out[1]:13.5f}{disp_out[2]:13.5f}\n")
        
        # Print moments
        fp.write("%\n%          CUMULATIVE MOMENTS\n")
        fp.write("%               MEMBER       END MOMENTS            NODES     PLASTIC MOM\n")
        
        for I in range(1, NM+1):
            k = 2*I - 1
            fp.write(f"%{I:19d}{CM[k]:14.2f}{CM[k+1]:11.2f}{MCON[I,1]:7d} AND{MCON[I,2]:2d}{PM[I]:14.2f}\n")
        
        # Print axial forces
        fp.write("%\n%          CUMULATIVE TENSION FORCES