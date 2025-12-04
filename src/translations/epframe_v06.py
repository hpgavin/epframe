#!/usr/bin/env python3
"""
ELASTIC PLASTIC ANALYSIS OF A PLANE FRAME
Translated from FORTRAN epframe.f by Hacksoo Lee, 1986
Python version with index-1 arrays using numpy

Nomenclature:
    Node (N) - connection point in the structure
    Element (E) - structural member connecting two nodes
"""

import numpy as np
import sys
from math import sqrt

class Index1Array:
    """Wrapper to provide 1-based indexing for numpy arrays"""
    def __init__(self, shape, dtype=float):
        if isinstance(shape, (int, np.integer)):
            self._data = np.zeros(shape + 1, dtype=dtype)
            self._shape = (shape,)
        elif isinstance(shape, (tuple, list)):
            if len(shape) == 1:
                self._data = np.zeros(shape[0] + 1, dtype=dtype)
                self._shape = (shape[0],)
            elif len(shape) == 2:
                self._data = np.zeros((shape[0] + 1, shape[1] + 1), dtype=dtype)
                self._shape = (shape[0], shape[1])
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
    
    def as_np(self):
        """Return view of data as 0-indexed numpy array (excludes index 0)"""
        if len(self._shape) == 1:
            return self._data[1:]
        else:
            return self._data[1:, 1:]
    
    def from_np(self, arr):
        """Copy data from 0-indexed numpy array into indices 1:n+1"""
        if len(self._shape) == 1:
            self._data[1:] = arr
        else:
            self._data[1:, 1:] = arr

def write_compact_header(fp, NCT, NE):
    """Write header line identifying columns in compact data record"""
    fp.write("%")
    fp.write(f"{'NCYCL':>4s}{'EL':>4s}{'ND':>4s}")
    fp.write(f"{'CLG':>16s}")
    # CX values - 3 per node (X, Y, Rotation)
    for nd in range(1, NCT+1):
        fp.write(f"{f'CX{nd}_X':>16s}")
        fp.write(f"{f'CX{nd}_Y':>16s}")
        fp.write(f"{f'CX{nd}_R':>16s}")
    # CM values - 2 per element (moment at each end)
    for el in range(1, NE+1):
        fp.write(f"{f'CM{el}_1':>16s}")
        fp.write(f"{f'CM{el}_2':>16s}")
    # CT values - 1 per element (axial force)
    for el in range(1, NE+1):
        fp.write(f"{f'CT{el}':>16s}")
    fp.write("\n")

def read_input_file(filename):
    """Read input file and return problem data"""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    idx = 0
    FN = int(lines[idx].split()[0])  # Frame number
    idx += 1
    
    parts = lines[idx].split()
    NCT, NE, E = int(parts[0]), int(parts[1]), float(parts[2])  # Node count, Number of elements, Modulus
    idx += 1
    
    # Node data
    CORD = Index1Array((NCT, 2))
    NTYPE = Index1Array((NCT, 3), dtype=int)
    
    for i in range(1, NCT+1):
        parts = lines[idx].split()
        node_num = int(parts[0])
        CORD[i,1] = float(parts[1])
        CORD[i,2] = float(parts[2])
        NTYPE[i,1] = int(parts[3])
        NTYPE[i,2] = int(parts[4])
        NTYPE[i,3] = int(parts[5])
        idx += 1
    
    # Element data
    ECON = Index1Array((NE, 2), dtype=int)  # Element connectivity
    SMA = Index1Array(NE)   # Second moment of area (I)
    AREA = Index1Array(NE)  # Cross-sectional area
    PM = Index1Array(NE)    # Plastic moment
    
    for i in range(1, NE+1):
        parts = lines[idx].split()
        elem_num = int(parts[0])
        ECON[i,1] = int(parts[1])
        ECON[i,2] = int(parts[2])
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
        ND = int(parts[0])  # Node number
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
    
    # Open output file
    fp = open(output_file, 'w')
    
    fp.write(f"%\n")
    fp.write(f"%     ELASTIC PLASTIC ANALYSIS OF FRAME NO {FN}\n")
    fp.write(f"%     ---------------------------------------\n%\n")
    
    # Calculate number of DOFs using numpy sum
    L = int(np.sum(NTYPE.as_np()))
    
    # Initialize load vector (zeros by default from Index1Array)
    VL = Index1Array(L)
    
    # Process loads
    OLEN = Index1Array(3)
    for load in loads:
        ND, FX, FY, FZ = load
        OLEN[1], OLEN[2], OLEN[3] = FX, FY, FZ
        
        LL = 0
        LN = ND - 1
        if LN > 0:
            for nd in range(1, LN+1):
                for k in range(1, 4):
                    LL += NTYPE[nd,k]
        
        for k in range(1, 4):
            if NTYPE[ND,k]:
                LL += 1
                VL[LL] = OLEN[k]
    
    fp.write(f"%     * GENERAL DATA\n")
    fp.write(f"%          NUMBER OF NODES         {NCT:6d}\n")
    fp.write(f"%          NUMBER OF ELEMENTS      {NE:6d}\n")
    fp.write(f"%          MOD OF ELASTICITY  {E:12.1f}\n%\n")
    
    # Input Echo: Node Data
    fp.write("%\n%     * DATA FOR NODES\n")
    fp.write("%           NODE   X-COORD   Y-COORD   DFX   DFY   DFZ\n%\n")
    for i in range(1, NCT+1):
        fp.write(f"%{i:14d}{CORD[i,1]:12.2f}{CORD[i,2]:10.2f}{NTYPE[i,1]:6d}{NTYPE[i,2]:6d}{NTYPE[i,3]:6d}\n")
    
    # Input Echo: Element Data
    fp.write("%\n%     * DATA FOR ELEMENTS\n")
    fp.write("%         ELEMENT    N1      N2       IXX      AREA        MP\n%\n")
    for i in range(1, NE+1):
        fp.write(f"%{i:14d}{ECON[i,1]:9d}{ECON[i,2]:8d}{SMA[i]:10.2f}{AREA[i]:10.2f}{PM[i]:10.2f}\n")
    
    # Input Echo: Load Data
    fp.write("%\n%     * DATA FOR LOADS\n")
    fp.write("%           NODE        PX        PY        PZ\n")
    for load in loads:
        ND, FX, FY, FZ = load
        fp.write(f"%{ND:14d}{FX:12.2f}{FY:10.2f}{FZ:10.2f}\n")
    
    # Calculate element lengths
    OLEN_elems = Index1Array(NE)
    for i in range(1, NE+1):
        N1 = ECON[i,1]
        N2 = ECON[i,2]
        X = CORD[N1,1] - CORD[N2,1]
        Y = CORD[N1,2] - CORD[N2,2]
        OLEN_elems[i] = sqrt(X*X + Y*Y)
    
    # Calculate stiffness coefficients
    SF = Index1Array((2*NE, 2))
    SA = Index1Array(NE)
    
    for i in range(1, NE+1):
        SF[2*i, 2] = 4.0 * E * SMA[i] / OLEN_elems[i]
        SF[2*i-1, 1] = SF[2*i, 2]
        SF[2*i, 1] = 0.5 * SF[2*i, 2]
        SF[2*i-1, 2] = SF[2*i, 1]
        SA[i] = E * AREA[i] / OLEN_elems[i]
    
    E2 = 2 * NE  # Number of element end moments
    E3 = 3 * NE  # Total element DOFs (2 moments + 1 axial per element)
    
    # Initialize cumulative variables (zeros by default from Index1Array)
    CM = Index1Array(E2)   # Cumulative moments
    CT = Index1Array(NE)   # Cumulative tensions
    CX = Index1Array(L)    # Cumulative displacements
    
    NCYCL = 0
    CLG = 0.0
    
    # Write initial state record (cycle 0)
    # Format: NCYCL, element, node, CLG, then CX values, CM values, CT values
    fp.write("%\n")
    write_compact_header(fp, NCT, NE)
    fp.write(f"{0:4d}{0:4d}{0:4d}")
    fp.write(f"{0.0:16.4E}")  # CLG = 0
    # CX values (all zeros) - need 3*NCT slots but only L are active
    for i in range(1, 3*NCT+1):
        fp.write(f"{0.0:16.4E}")
    # CM values (2*NE)
    for i in range(1, E2+1):
        fp.write(f"{CM[i]:16.4E}")
    # CT values (NE)
    for i in range(1, NE+1):
        fp.write(f"{CT[i]:16.4E}")
    fp.write("\n")
    
    # Build flexibility matrix K (compatibility matrix) - zeros by default
    K = Index1Array((L, E3))
    
    row_count = 0
    row_max = 0
    
    for nd in range(1, NCT+1):
        for el in range(1, NE+1):
            NA = row_count
            
            # Check if node nd is connected to element el
            if nd == ECON[el,1]:
                NF = ECON[el,2]  # Far node
                EI_near = 2*el - 1  # Element moment index at near end
                EI_far = EI_near + 1
            elif nd == ECON[el,2]:
                NF = ECON[el,1]
                EI_near = 2*el
                EI_far = EI_near - 1
            else:
                continue
            
            X = CORD[NF,1] - CORD[nd,1]
            Y = CORD[NF,2] - CORD[nd,2]
            D = sqrt(X*X + Y*Y)
            S = Y / D
            C = X / D
            axial_col = 2*NE + el  # Column for axial force
            
            if NTYPE[nd,1]:
                NA = NA + 1
                K[NA,EI_near] = S / D
                K[NA,EI_far] = K[NA,EI_near]
                K[NA,axial_col] = -C
            
            if NTYPE[nd,2]:
                NA = NA + 1
                K[NA,EI_near] = -C / D
                K[NA,EI_far] = K[NA,EI_near]
                K[NA,axial_col] = -S
            
            if NTYPE[nd,3]:
                NA = NA + 1
                K[NA,EI_near] = 1.0
            
            if NA > row_max:
                row_max = NA
        
        row_count = row_max
    
    # Main analysis loop
    while True:
        NCYCL += 1
        
        # Build the full E3 x E3 element stiffness matrix S_full
        # Structure: [SF_block (E2xE2) | 0      ]
        #            [0                | SA_diag]
        S_full = np.zeros((E3, E3))
        
        # Fill in 2x2 blocks for each element's flexural stiffness
        for i in range(1, NE+1):
            row = 2*i - 2  # 0-indexed
            S_full[row, row] = SF[2*i-1, 1]
            S_full[row, row+1] = SF[2*i-1, 2]
            S_full[row+1, row] = SF[2*i, 1]
            S_full[row+1, row+1] = SF[2*i, 2]
        
        # Fill diagonal for axial stiffness
        for i in range(1, NE+1):
            S_full[E2-1+i, E2-1+i] = SA[i]
        
        # Get 0-indexed numpy views
        K_np = K.as_np()  # L x E3
        VL_np = VL.as_np()  # L
        
        # Form stiffness matrix: KSAT = K @ S_full @ K.T
        KSAT_np = K_np @ S_full @ K_np.T
        
        # Solve system using numpy linear algebra
        try:
            disp_np = np.linalg.solve(KSAT_np, VL_np)
        except np.linalg.LinAlgError:
            fp.write("%     * DIVISION BY ZERO IN SOLUTION OF EQUATION\n%\n")
            break
        
        # Check for excessive deformations
        XLMT = 1000.0
        max_disp = np.max(np.abs(disp_np))
        
        if max_disp > XLMT:
            fp.write(f"%\n%     *** DEFORMATIONS LARGER THAN {XLMT:.1f} IN CYCLE NO {NCYCL:4d}\n%\n")
            break
        
        # Calculate element deformations: CSAT = K.T @ disp
        CSAT_np = K_np.T @ disp_np
        
        # Calculate element forces: SATX = S_full @ CSAT
        SATX_np = S_full @ CSAT_np
        
        # Convert back to 1-indexed for remaining calculations
        disp = Index1Array(L)
        disp.from_np(disp_np)
        CSAT = Index1Array(E3)
        CSAT.from_np(CSAT_np)
        SATX = Index1Array(E3)
        SATX.from_np(SATX_np)
        
        # Calculate load factors to plastic hinge using numpy
        # Get numpy views
        SATX_e2 = SATX.as_np()[:E2]  # First E2 values (moments)
        CM_np = CM.as_np()
        PM_np = PM.as_np()
        
        # For each moment location I, the element is k = (I+1)//2 (1-indexed)
        # In 0-indexed: elem_idx = I // 2 for I in 0..E2-1
        elem_indices = np.arange(E2) // 2  # 0-indexed element for each moment location
        PM_expanded = PM_np[elem_indices]  # Plastic moment for each location
        
        # Calculate ZERO threshold and ALG
        ZERO = 0.001 * PM_expanded
        small_satx = np.abs(SATX_e2) < ZERO
        # Use safe division: replace small values with 1 to avoid divide-by-zero
        safe_denom = np.where(small_satx, 1.0, np.abs(SATX_e2))
        ALG_np = np.where(small_satx, 
                         1.0E10, 
                         (PM_expanded - np.abs(CM_np)) / safe_denom)
        
        # Find minimum load factor where TEST >= 0
        TEST = CM_np * SATX_e2
        valid_mask = TEST >= 0.0
        ALG_masked = np.where(valid_mask, ALG_np, 1.0E10)
        NPH = int(np.argmin(ALG_masked)) + 1  # Convert to 1-indexed
        SALG = ALG_masked[NPH - 1]
        
        # Scale forces and displacements, update cumulative values using numpy
        SATX_np *= SALG
        SATX.from_np(SATX_np)
        
        CLG += SALG
        
        # Update CM (first E2 of SATX)
        CM_np += SATX_np[:E2]
        CM.from_np(CM_np)
        
        # Update CT (last NE of SATX)
        CT.as_np()[:] += SATX_np[E2:]
        
        # Update CX
        disp_np *= SALG
        CX.as_np()[:] += disp_np
        
        # Determine which element and node
        EL = (NPH + 1) // 2
        k = (NPH // 2) * 2 - NPH
        if k:
            ND = ECON[EL,1]
        else:
            ND = ECON[EL,2]
        
        EL_hinge = EL
        ND_hinge = ND
        
        # Print results
        fp.write("%\n%\n%\n")
        fp.write(f"%     * PLASTIC HINGE {NCYCL:3d} FORMED IN ELEMENT {EL:3d} NEAR NODE {ND:3d} WHEN LOAD FACTOR IS {CLG:12.3f}\n%\n")
        
        fp.write("%          CUMULATIVE DEFORMATIONS\n")
        fp.write("%                NODE    X-DISP       Y-DISP       ROTN\n")
        
        LL = 0
        for nd in range(1, NCT+1):
            disp_out = [0.0, 0.0, 0.0]
            for dof in range(1, 4):
                if NTYPE[nd,dof]:
                    LL += 1
                    disp_out[dof-1] = CX[LL]
            fp.write(f"%{nd:19d}{disp_out[0]:13.5f}{disp_out[1]:13.5f}{disp_out[2]:13.5f}\n")
        
        fp.write("%\n%          CUMULATIVE MOMENTS\n")
        fp.write("%             ELEMENT       END MOMENTS             NODES     PLASTIC MOM\n")
        
        for el in range(1, NE+1):
            k = 2*el - 1
            fp.write(f"%{el:19d}{CM[k]:14.2f}{CM[k+1]:11.2f}{ECON[el,1]:7d} AND{ECON[el,2]:2d}{PM[el]:14.2f}\n")
        
        fp.write("%\n%          CUMULATIVE TENSION FORCES\n")
        fp.write("%             ELEMENT     TENSION\n")
        
        for el in range(1, NE+1):
            fp.write(f"%{el:19d}{CT[el]:15.2f}\n")
        
        # Calculate and output reactions at support nodes
        # Reactions are what the support provides (opposite of element forces on node)
        fp.write("%\n%          REACTIONS AT SUPPORTS\n")
        fp.write("%                NODE       RX           RY           MZ\n")
        
        for nd in range(1, NCT+1):
            # Check if this node has any restrained DOFs (is a support)
            is_support = (NTYPE[nd,1] == 0 or NTYPE[nd,2] == 0 or NTYPE[nd,3] == 0)
            if not is_support:
                continue
            
            Rx, Ry, Rz = 0.0, 0.0, 0.0
            
            # Sum contributions from all elements connected to this node
            for el in range(1, NE+1):
                N1 = ECON[el,1]
                N2 = ECON[el,2]
                
                if N1 != nd and N2 != nd:
                    continue
                
                # Get element geometry (direction from N1 to N2)
                dx = CORD[N2,1] - CORD[N1,1]
                dy = CORD[N2,2] - CORD[N1,2]
                L_el = sqrt(dx*dx + dy*dy)
                Cx = dx / L_el  # cosine
                Sy = dy / L_el  # sine
                
                # Get element forces
                Me1 = CM[2*el - 1]  # Moment at N1 end
                Me2 = CM[2*el]      # Moment at N2 end
                N = CT[el]          # Axial force (tension positive)
                V = (Me1 + Me2) / L_el  # Shear (from moment equilibrium)
                
                if N1 == nd:
                    # Element forces on node at end 1
                    Rx += N * Cx - V * Sy
                    Ry += N * Sy + V * Cx
                    Rz += Me1
                else:  # N2 == nd
                    # Element forces on node at end 2
                    Rx += -N * Cx + V * Sy
                    Ry += -N * Sy - V * Cx
                    Rz += Me2
            
            # Reaction is negative of element forces (what support provides)
            fp.write(f"%{nd:19d}{-Rx:13.2f}{-Ry:13.2f}{-Rz:13.2f}\n")
        
        # Write compact data record for post-processing
        # Format: NCYCL, element, node, CLG, then CX values, CM values, CT values
        fp.write("%\n")
        write_compact_header(fp, NCT, NE)
        fp.write(f"{NCYCL:4d}{EL_hinge:4d}{ND_hinge:4d}")
        fp.write(f"{CLG:16.4E}")
        # CX values - write for all 3*NCT slots (zeros for restrained DOFs)
        LL = 0
        for nd in range(1, NCT+1):
            for dof in range(1, 4):
                if NTYPE[nd,dof]:
                    LL += 1
                    fp.write(f"{CX[LL]:16.4E}")
                else:
                    fp.write(f"{0.0:16.4E}")
        # CM values (2*NE)
        for i in range(1, E2+1):
            fp.write(f"{CM[i]:16.4E}")
        # CT values (NE)
        for i in range(1, NE+1):
            fp.write(f"{CT[i]:16.4E}")
        fp.write("\n")
        
        # Modify stiffness for plastic hinge
        ITEST = ((NPH//2)*2) - NPH
        if ITEST:
            SF[NPH+1, 2] = 0.75 * SF[NPH+1, 2]
            SF[NPH+1, 1] = 0.0
            SF[NPH, 1] = 0.0
            SF[NPH, 2] = 0.0
        else:
            SF[NPH-1, 1] = 0.75 * SF[NPH-1, 1]
            SF[NPH-1, 2] = 0.0
            SF[NPH, 1] = 0.0
            SF[NPH, 2] = 0.0
    
    fp.write(f"%\n%     ANALYSIS COMPLETED FOR FRAME NO {FN:3d}\n\n")
    fp.close()
    print(f"Analysis completed. Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python epframe.py input_file output_file")
        sys.exit(1)
    
    epframe_analysis(sys.argv[1], sys.argv[2])
