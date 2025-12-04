#!/usr/bin/env python3
"""
ELASTIC PLASTIC ANALYSIS OF A PLANE FRAME
Translated from FORTRAN epframe.f by Hacksoo Lee, 1986
Python version using standard 0-indexed numpy arrays

Nomenclature:
    Node (N) - connection point in the structure (0-indexed internally, 1-indexed in I/O)
    Element (E) - structural member connecting two nodes (0-indexed internally, 1-indexed in I/O)
"""

import numpy as np
import sys
from math import sqrt

def get_csv_header(NCT, NE):
    """Generate CSV header row"""
    headers = ['NCYCL', 'EL', 'ND', 'CLG']
    # CX values - 3 per node (X, Y, Rotation)
    for nd in range(NCT):
        headers.append(f'CX{nd+1}_X')
        headers.append(f'CX{nd+1}_Y')
        headers.append(f'CX{nd+1}_R')
    # CM values - 2 per element (moment at each end)
    for el in range(NE):
        headers.append(f'CM{el+1}_1')
        headers.append(f'CM{el+1}_2')
    # CT values - 1 per element (axial force)
    for el in range(NE):
        headers.append(f'CT{el+1}')
    return ','.join(headers)

def write_csv_row(csv_fp, NCYCL, EL, ND, CLG, CX, CM, CT, NTYPE, NCT, NE):
    """Write a data row to CSV file"""
    values = [f'{NCYCL}', f'{EL}', f'{ND}', f'{CLG:.6E}']
    # CX values - write for all 3*NCT slots
    LL = 0
    for nd in range(NCT):
        for dof in range(3):
            if NTYPE[nd, dof]:
                values.append(f'{CX[LL]:.6E}')
                LL += 1
            else:
                values.append(f'{0.0:.6E}')
    # CM values
    for i in range(2*NE):
        values.append(f'{CM[i]:.6E}')
    # CT values
    for i in range(NE):
        values.append(f'{CT[i]:.6E}')
    csv_fp.write(','.join(values) + '\n')

def read_input_file(filename):
    """Read input file and return problem data
    
    Comments: Lines starting with # or text after # are ignored
    All indices are converted to 0-based internally
    """
    with open(filename, 'r') as f:
        lines = []
        for line in f:
            # Remove comments (everything from # to end of line)
            line = line.split('#')[0].strip()
            if line:  # Skip empty lines
                lines.append(line)
    
    idx = 0
    FN = int(lines[idx].split()[0])  # Frame number
    idx += 1
    
    parts = lines[idx].split()
    NCT, NE, E = int(parts[0]), int(parts[1]), float(parts[2])  # Node count, Number of elements, Modulus
    idx += 1
    
    # Node data (0-indexed)
    CORD = np.zeros((NCT, 2))
    NTYPE = np.zeros((NCT, 3), dtype=int)
    
    for i in range(NCT):
        parts = lines[idx].split()
        # Input node number is 1-based, but we store at 0-based index
        node_num = int(parts[0]) - 1  # Convert to 0-indexed
        CORD[node_num, 0] = float(parts[1])
        CORD[node_num, 1] = float(parts[2])
        NTYPE[node_num, 0] = int(parts[3])
        NTYPE[node_num, 1] = int(parts[4])
        NTYPE[node_num, 2] = int(parts[5])
        idx += 1
    
    # Element data (0-indexed, connectivity stores 0-based node indices)
    ECON = np.zeros((NE, 2), dtype=int)  # Element connectivity
    SMA = np.zeros(NE)   # Second moment of area (I)
    AREA = np.zeros(NE)  # Cross-sectional area
    PM = np.zeros(NE)    # Plastic moment
    
    for i in range(NE):
        parts = lines[idx].split()
        elem_num = int(parts[0]) - 1  # Convert to 0-indexed
        ECON[elem_num, 0] = int(parts[1]) - 1  # Node 1 (0-indexed)
        ECON[elem_num, 1] = int(parts[2]) - 1  # Node 2 (0-indexed)
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
        ND = int(parts[0]) - 1  # Convert to 0-indexed
        FX = float(parts[1])
        FY = float(parts[2])
        FZ = float(parts[3])
        loads.append((ND, FX, FY, FZ))
        idx += 1
    
    return FN, NCT, NE, E, CORD, NTYPE, ECON, SMA, AREA, PM, loads

def epframe_analysis(input_file, output_file):
    """Main elastic-plastic frame analysis"""
    
    # Read input
    FN, NCT, NE, E, CORD, NTYPE, ECON, SMA, AREA, PM, loads = read_input_file(input_file)
    
    # Open output file and CSV file
    fp = open(output_file, 'w')
    csv_file = output_file + '.csv'
    csv_fp = open(csv_file, 'w')
    
    # Write CSV header
    csv_fp.write(get_csv_header(NCT, NE) + '\n')
    
    fp.write(f"%\n")
    fp.write(f"%     ELASTIC PLASTIC ANALYSIS OF FRAME NO {FN}\n")
    fp.write(f"%     ---------------------------------------\n%\n")
    
    # Calculate number of DOFs
    L = int(np.sum(NTYPE))
    
    # Initialize load vector
    VL = np.zeros(L)
    
    # Process loads
    for load in loads:
        ND, FX, FY, FZ = load
        OLEN = [FX, FY, FZ]
        
        # Count DOFs before this node
        LL = 0
        for nd in range(ND):
            LL += np.sum(NTYPE[nd])
        
        # Assign loads to corresponding DOFs
        for k in range(3):
            if NTYPE[ND, k]:
                VL[LL] = OLEN[k]
                LL += 1
    
    fp.write(f"%     * GENERAL DATA\n")
    fp.write(f"%          NUMBER OF NODES         {NCT:6d}\n")
    fp.write(f"%          NUMBER OF ELEMENTS      {NE:6d}\n")
    fp.write(f"%          MOD OF ELASTICITY  {E:12.1f}\n%\n")
    
    # Input Echo: Node Data
    fp.write("%\n%     * DATA FOR NODES\n")
    fp.write("%           NODE   X-COORD   Y-COORD   DFX   DFY   DFZ\n%\n")
    for i in range(NCT):
        fp.write(f"%{i+1:14d}{CORD[i,0]:12.2f}{CORD[i,1]:10.2f}{NTYPE[i,0]:6d}{NTYPE[i,1]:6d}{NTYPE[i,2]:6d}\n")
    
    # Input Echo: Element Data
    fp.write("%\n%     * DATA FOR ELEMENTS\n")
    fp.write("%         ELEMENT    N1      N2       IXX      AREA        MP\n%\n")
    for i in range(NE):
        fp.write(f"%{i+1:14d}{ECON[i,0]+1:9d}{ECON[i,1]+1:8d}{SMA[i]:10.2f}{AREA[i]:10.2f}{PM[i]:10.2f}\n")
    
    # Input Echo: Load Data
    fp.write("%\n%     * DATA FOR LOADS\n")
    fp.write("%           NODE        PX        PY        PZ\n")
    for load in loads:
        ND, FX, FY, FZ = load
        fp.write(f"%{ND+1:14d}{FX:12.2f}{FY:10.2f}{FZ:10.2f}\n")
    
    # Calculate element lengths
    OLEN_elems = np.zeros(NE)
    for i in range(NE):
        N1, N2 = ECON[i, 0], ECON[i, 1]
        dx = CORD[N1, 0] - CORD[N2, 0]
        dy = CORD[N1, 1] - CORD[N2, 1]
        OLEN_elems[i] = sqrt(dx*dx + dy*dy)
    
    # Calculate stiffness coefficients
    # SF[2*i, :] and SF[2*i+1, :] are the 2x2 flexural stiffness for element i
    SF = np.zeros((2*NE, 2))
    SA = np.zeros(NE)
    
    for i in range(NE):
        SF[2*i+1, 1] = 4.0 * E * SMA[i] / OLEN_elems[i]
        SF[2*i, 0] = SF[2*i+1, 1]
        SF[2*i+1, 0] = 0.5 * SF[2*i+1, 1]
        SF[2*i, 1] = SF[2*i+1, 0]
        SA[i] = E * AREA[i] / OLEN_elems[i]
    
    E2 = 2 * NE  # Number of element end moments
    E3 = 3 * NE  # Total element DOFs (2 moments + 1 axial per element)
    
    # Initialize cumulative variables
    CM = np.zeros(E2)   # Cumulative moments
    CT = np.zeros(NE)   # Cumulative tensions
    CX = np.zeros(L)    # Cumulative displacements
    
    NCYCL = 0
    CLG = 0.0
    
    # Write initial state record (cycle 0) to CSV
    write_csv_row(csv_fp, 0, 0, 0, 0.0, CX, CM, CT, NTYPE, NCT, NE)
    
    # Build compatibility matrix K
    # Each row corresponds to a DOF, each column to an element force
    K = np.zeros((L, E3))
    
    row_start = 0  # Starting row for current node
    
    for nd in range(NCT):
        for el in range(NE):
            # Check if node nd is connected to element el
            if nd == ECON[el, 0]:
                NF = ECON[el, 1]  # Far node
                EI_near = 2*el      # Element moment index at near end (0-indexed)
                EI_far = 2*el + 1
            elif nd == ECON[el, 1]:
                NF = ECON[el, 0]
                EI_near = 2*el + 1
                EI_far = 2*el
            else:
                continue
            
            X = CORD[NF, 0] - CORD[nd, 0]
            Y = CORD[NF, 1] - CORD[nd, 1]
            D = sqrt(X*X + Y*Y)
            S = Y / D
            C = X / D
            axial_col = 2*NE + el  # Column for axial force
            
            # Fill rows for this node's free DOFs
            NA = row_start
            if NTYPE[nd, 0]:
                K[NA, EI_near] = S / D
                K[NA, EI_far] = K[NA, EI_near]
                K[NA, axial_col] = -C
                NA += 1
            
            if NTYPE[nd, 1]:
                K[NA, EI_near] = -C / D
                K[NA, EI_far] = K[NA, EI_near]
                K[NA, axial_col] = -S
                NA += 1
            
            if NTYPE[nd, 2]:
                K[NA, EI_near] = 1.0
        
        # Advance row_start by number of DOFs at this node
        row_start += int(np.sum(NTYPE[nd]))
    
    # Main analysis loop
    while True:
        NCYCL += 1
        
        # Build the full E3 x E3 element stiffness matrix S_full
        S_full = np.zeros((E3, E3))
        
        # Fill in 2x2 blocks for each element's flexural stiffness
        for i in range(NE):
            row = 2*i
            S_full[row, row] = SF[2*i, 0]
            S_full[row, row+1] = SF[2*i, 1]
            S_full[row+1, row] = SF[2*i+1, 0]
            S_full[row+1, row+1] = SF[2*i+1, 1]
        
        # Fill diagonal for axial stiffness
        for i in range(NE):
            S_full[E2+i, E2+i] = SA[i]
        
        # Form stiffness matrix: KSAT = K @ S_full @ K.T
        KSAT = K @ S_full @ K.T
        
        # Solve system using numpy linear algebra
        try:
            disp = np.linalg.solve(KSAT, VL)
        except np.linalg.LinAlgError:
            fp.write("%     * DIVISION BY ZERO IN SOLUTION OF EQUATION\n%\n")
            break
        
        # Check for excessive deformations
        XLMT = 1000.0
        max_disp = np.max(np.abs(disp))
        
        if max_disp > XLMT:
            fp.write(f"%\n%     *** DEFORMATIONS LARGER THAN {XLMT:.1f} IN CYCLE NO {NCYCL:4d}\n%\n")
            break
        
        # Calculate element deformations and forces
        CSAT = K.T @ disp
        SATX = S_full @ CSAT
        
        # Calculate load factors to plastic hinge
        SATX_e2 = SATX[:E2]  # Moments
        
        # Element index for each moment location
        elem_indices = np.arange(E2) // 2
        PM_expanded = PM[elem_indices]
        
        # Calculate load factor to reach plastic moment
        ZERO = 0.001 * PM_expanded
        small_satx = np.abs(SATX_e2) < ZERO
        safe_denom = np.where(small_satx, 1.0, np.abs(SATX_e2))
        ALG = np.where(small_satx, 1.0E10, (PM_expanded - np.abs(CM)) / safe_denom)
        
        # Find minimum load factor where moment is increasing (same sign)
        TEST = CM * SATX_e2
        valid_mask = TEST >= 0.0
        ALG_masked = np.where(valid_mask, ALG, 1.0E10)
        NPH = int(np.argmin(ALG_masked))  # 0-indexed
        SALG = ALG_masked[NPH]
        
        # Scale forces and update cumulative values
        SATX *= SALG
        CLG += SALG
        CM += SATX[:E2]
        CT += SATX[E2:]
        disp *= SALG
        CX += disp
        
        # Determine which element and node
        EL = NPH // 2  # 0-indexed element
        if NPH % 2 == 0:
            ND = ECON[EL, 0]
        else:
            ND = ECON[EL, 1]
        
        EL_hinge = EL
        ND_hinge = ND
        
        # Print results (convert to 1-indexed for display)
        fp.write("%\n%\n%\n")
        fp.write(f"%     * PLASTIC HINGE {NCYCL:3d} FORMED IN ELEMENT {EL+1:3d} NEAR NODE {ND+1:3d} WHEN LOAD FACTOR IS {CLG:12.3f}\n%\n")
        
        fp.write("%          CUMULATIVE DEFORMATIONS\n")
        fp.write("%                NODE    X-DISP       Y-DISP       ROTN\n")
        
        LL = 0
        for nd in range(NCT):
            disp_out = [0.0, 0.0, 0.0]
            for dof in range(3):
                if NTYPE[nd, dof]:
                    disp_out[dof] = CX[LL]
                    LL += 1
            fp.write(f"%{nd+1:19d}{disp_out[0]:13.5f}{disp_out[1]:13.5f}{disp_out[2]:13.5f}\n")
        
        fp.write("%\n%          CUMULATIVE MOMENTS\n")
        fp.write("%             ELEMENT       END MOMENTS             NODES     PLASTIC MOM\n")
        
        for el in range(NE):
            k = 2*el
            fp.write(f"%{el+1:19d}{CM[k]:14.2f}{CM[k+1]:11.2f}{ECON[el,0]+1:7d} AND{ECON[el,1]+1:2d}{PM[el]:14.2f}\n")
        
        fp.write("%\n%          CUMULATIVE TENSION FORCES\n")
        fp.write("%             ELEMENT     TENSION\n")
        
        for el in range(NE):
            fp.write(f"%{el+1:19d}{CT[el]:15.2f}\n")
        
        # Calculate and output reactions at support nodes
        fp.write("%\n%          REACTIONS AT SUPPORTS\n")
        fp.write("%                NODE       RX           RY           MZ\n")
        
        for nd in range(NCT):
            # Check if this node has any restrained DOFs
            is_support = (NTYPE[nd, 0] == 0 or NTYPE[nd, 1] == 0 or NTYPE[nd, 2] == 0)
            if not is_support:
                continue
            
            Rx, Ry, Rz = 0.0, 0.0, 0.0
            
            # Sum contributions from all elements connected to this node
            for el in range(NE):
                N1, N2 = ECON[el, 0], ECON[el, 1]
                
                if N1 != nd and N2 != nd:
                    continue
                
                # Get element geometry (direction from N1 to N2)
                dx = CORD[N2, 0] - CORD[N1, 0]
                dy = CORD[N2, 1] - CORD[N1, 1]
                L_el = sqrt(dx*dx + dy*dy)
                Cx = dx / L_el
                Sy = dy / L_el
                
                # Get element forces
                Me1 = CM[2*el]      # Moment at N1 end
                Me2 = CM[2*el + 1]  # Moment at N2 end
                N_ax = CT[el]       # Axial force (tension positive)
                V = (Me1 + Me2) / L_el  # Shear
                
                if N1 == nd:
                    Rx += N_ax * Cx - V * Sy
                    Ry += N_ax * Sy + V * Cx
                    Rz += Me1
                else:
                    Rx += -N_ax * Cx + V * Sy
                    Ry += -N_ax * Sy - V * Cx
                    Rz += Me2
            
            fp.write(f"%{nd+1:19d}{-Rx:13.2f}{-Ry:13.2f}{-Rz:13.2f}\n")
        
        # Write data row to CSV
        write_csv_row(csv_fp, NCYCL, EL_hinge+1, ND_hinge+1, CLG, CX, CM, CT, NTYPE, NCT, NE)
        
        # Modify stiffness for plastic hinge
        if NPH % 2 == 0:
            # Hinge at first end of element
            SF[NPH+1, 1] = 0.75 * SF[NPH+1, 1]
            SF[NPH+1, 0] = 0.0
            SF[NPH, 0] = 0.0
            SF[NPH, 1] = 0.0
        else:
            # Hinge at second end of element
            SF[NPH-1, 0] = 0.75 * SF[NPH-1, 0]
            SF[NPH-1, 1] = 0.0
            SF[NPH, 0] = 0.0
            SF[NPH, 1] = 0.0
    
    fp.write(f"%\n%     ANALYSIS COMPLETED FOR FRAME NO {FN:3d}\n\n")
    fp.close()
    csv_fp.close()
    print(f"Analysis completed. Results written to {output_file}")
    print(f"Compact data written to {csv_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python epframe.py input_file output_file")
        sys.exit(1)
    
    epframe_analysis(sys.argv[1], sys.argv[2])
