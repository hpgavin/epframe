#!/usr/bin/env python3
"""
ELASTIC PLASTIC ANALYSIS OF A PLANE FRAME
Translated from FORTRAN epframe.f by Hacksoo Lee, 1986
Python version with index-1 arrays using numpy
"""

import numpy as np
import sys
from math import sqrt

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

def lu_decompose(A, n):
    """
    LU decomposition with partial pivoting - optimized with numpy
    Returns L, U, and permutation vector
    """
    # Extract numpy array from Index1Array
    A_np = A._data[1:n+1, 1:n+1].copy()
    
    perm = np.arange(n)
    
    for k in range(n):
        # Find pivot using numpy
        pivot_row = k + np.argmax(np.abs(A_np[k:, k]))
        
        # Swap rows if necessary
        if pivot_row != k:
            A_np[[k, pivot_row]] = A_np[[pivot_row, k]]
            perm[[k, pivot_row]] = perm[[pivot_row, k]]
        
        # Check for singular matrix
        if abs(A_np[k, k]) < 1e-20:
            return None, None, None
        
        # Elimination using vectorized operations
        if k < n - 1:
            A_np[k+1:, k] /= A_np[k, k]
            A_np[k+1:, k+1:] -= np.outer(A_np[k+1:, k], A_np[k, k+1:])
    
    # Create Index1Array results
    LU = Index1Array((n, n))
    perm_out = Index1Array(n, dtype=int)
    
    for i in range(n):
        perm_out[i+1] = perm[i] + 1  # Convert to 1-indexed
        for j in range(n):
            LU[i+1, j+1] = A_np[i, j]
    
    return LU, LU, perm_out  # Return same matrix for both L and U (combined storage)

def lu_solve(A, b, n):
    """
    Solve Ax = b using LU decomposition - optimized with numpy
    Returns solution vector x (1-indexed)
    """
    LU, _, perm = lu_decompose(A, n)
    if LU is None:
        return None
    
    # Extract numpy arrays
    LU_np = LU._data[1:n+1, 1:n+1]
    b_np = b._data[1:n+1].copy()
    perm_np = perm._data[1:n+1] - 1  # Convert to 0-indexed
    
    # Permute b according to pivot order
    b_np = b_np[perm_np]
    
    # Forward substitution: Ly = Pb (vectorized)
    y = np.zeros(n)
    for i in range(n):
        y[i] = b_np[i] - LU_np[i, :i] @ y[:i]
    
    # Back substitution: Ux = y (vectorized with @)
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - LU_np[i, i+1:] @ x[i+1:]) / LU_np[i, i]
    
    # Convert to Index1Array
    x_out = Index1Array(n)
    for i in range(n):
        x_out[i+1] = x[i]
    
    return x_out

def parse_fixity(fixity_str):
    """
    Parse fixity string to DOF flags
    Input: string like 'X Y Z', '* Y *', 'x y z', etc.
    Returns: (DFX, DFY, DFZ) where 1=free, 0=fixed
    """
    fixity_str = fixity_str.upper().strip()
    parts = fixity_str.split()
    
    # Initialize all as free (1)
    dof = {'X': 1, 'Y': 1, 'Z': 1}
    
    # Mark as fixed (0) if letter appears
    for part in parts:
        if part in ['X', 'Y', 'Z']:
            dof[part] = 0
    
    return dof['X'], dof['Y'], dof['Z']

def read_input_file(filename):
    """Read input file with comment support and return problem data"""
    with open(filename, 'r') as f:
        # Remove comments (everything after #) and blank lines
        lines = []
        for line in f:
            # Remove comment
            line = line.split('#')[0].strip()
            if line:
                lines.append(line)
    
    idx = 0
    JFN = int(lines[idx].split()[0])
    idx += 1
    
    parts = lines[idx].split()
    NM, JCT, E = int(parts[0]), int(parts[1]), float(parts[2])
    idx += 1
    
    # Node data
    CORD = Index1Array((JCT, 2))
    JTYPE = Index1Array((JCT, 3), dtype=int)
    
    for i in range(1, JCT+1):
        parts = lines[idx].split()
        node_num = int(parts[0])
        CORD[i,1] = float(parts[1])
        CORD[i,2] = float(parts[2])
        
        # Check if old format (numeric) or new format (fixity string)
        if len(parts) >= 6 and parts[3].replace('.','').replace('-','').isdigit():
            # Old format: I X Y DFX DFY DFZ (numeric)
            JTYPE[i,1] = int(parts[3])
            JTYPE[i,2] = int(parts[4])
            JTYPE[i,3] = int(parts[5])
        else:
            # New format: N X Y fixity_string
            # Combine remaining parts as fixity string
            fixity_str = ' '.join(parts[3:])
            JTYPE[i,1], JTYPE[i,2], JTYPE[i,3] = parse_fixity(fixity_str)
        
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
    
    return JFN, JCT, NM, E, CORD, JTYPE, MCON, SMA, AREA, PM, loads

def epframe_analysis(input_file, output_file):
    """Main elastic-plastic frame analysis"""
    
    # Read input
    JFN, JCT, NM, E, CORD, JTYPE, MCON, SMA, AREA, PM, loads = read_input_file(input_file)
    
    # Open output file
    fp = open(output_file, 'w')
    
    fp.write(f"%\n")
    fp.write(f"%     ELASTIC PLASTIC ANALYSIS OF FRAME NO {JFN}\n")
    fp.write(f"%     ---------------------------------------\n%\n")
    
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
    
    fp.write(f"%     * GENERAL DATA\n")
    fp.write(f"%          NUMBER OF JOINTS        {JCT:6d}\n")
    fp.write(f"%          NUMBER OF MEMBERS       {NM:6d}\n")
    fp.write(f"%          MOD OF ELASTICITY  {E:12.1f}\n%\n")
    
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
    
    # Build flexibility matrix K (compatibility matrix)
    K = Index1Array((L, M3))
    for i in range(1, L+1):
        for j in range(1, M3+1):
            K[i,j] = 0.0
    
    NJ = 0
    NK = 0
    
    for J in range(1, JCT+1):
        for M in range(1, NM+1):
            NA = NJ
            
            # Check if joint J is connected to member M
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
    
    # Convert Index1Array to numpy for vectorized operations
    K_np = K._data[1:L+1, 1:M3+1]  # Extract actual data excluding index 0
    SF_np = SF._data[1:M2+1, 1:3]
    SA_np = SA._data[1:NM+1]
    
    # Main analysis loop
    while True:
        NCYCL += 1
        
        # Form stiffness matrix KSAT using vectorized operations
        # Build stiffness matrix S as a diagonal-like operation
        S_diag = np.zeros(M3)
        
        for I in range(M2):
            k = ((I+2)//2)*2 - 2  # Adjusted for 0-indexing
            S_diag[I] = SF_np[I, 0] if (I % 2 == 0) else SF_np[I, 0]
            if k+1 < M2:
                S_diag[I] += SF_np[I, 1] if (I % 2 == 1) else SF_np[I, 1]
        
        # Simpler approach: compute S*K row by row for flexural terms
        SK = np.zeros((L, M3))
        for j in range(L):
            for i in range(M2):
                k = ((i+1)//2)*2 - 1
                SK[j, i] = SF_np[i, 0]*K_np[j, k-1] + SF_np[i, 1]*K_np[j, k]
            for i in range(NM):
                SK[j, M2+i] = SA_np[i] * K_np[j, M2+i]
        
        # Form KSAT = K^T @ (S*K) using @ operator
        KSAT_np = K_np.T @ SK
        
        # Convert back to Index1Array format
        KSAT = Index1Array((L, L+1))
        for i in range(1, L+1):
            for j in range(1, L+1):
                KSAT[i,j] = KSAT_np[i-1, j-1]
        
        # Add load vector
        for I in range(1, L+1):
            KSAT[I,L+1] = VL[I]
        
        # Solve system using LU decomposition
        b = Index1Array(L)
        for i in range(1, L+1):
            b[i] = KSAT[i,L+1]
        
        disp = lu_solve(KSAT, b, L)
        
        if disp is None:
            fp.write("%     * DIVISION BY ZERO IN SOLUTION OF EQUATION\n%\n")
            break
        
        # Check for excessive deformations
        XLMT = 1000.0
        disp_np = disp._data[1:L+1]
        max_disp = np.max(np.abs(disp_np))
        
        if max_disp > XLMT:
            fp.write(f"%\n%     *** DEFORMATIONS LARGER THAN {XLMT:.1f} IN CYCLE NO {NCYCL:4d}\n%\n")
            break
        
        # Calculate member forces: CSAT = K^T @ disp using @ operator
        CSAT_np = K_np.T @ disp_np
        
        # Calculate moments and axial forces using vectorized operations
        SATX = Index1Array(M3)
        for I in range(1, M2+1):
            k = ((I+1)//2)*2 - 1
            SATX[I] = SF[I,1]*CSAT_np[k-1] + SF[I,2]*CSAT_np[k]
        
        for I in range(1, NM+1):
            k = M2 + I
            SATX[k] = SA[I] * CSAT_np[k-1]
        
        # Calculate load factors to plastic hinge - vectorized where possible
        ALG = Index1Array(M2)
        PM_np = PM._data[1:NM+1]
        CM_np = CM._data[1:M2+1]
        SATX_np = SATX._data[1:M2+1]
        
        for I in range(1, M2+1):
            k = (I+1)//2
            ZERO = 0.001 * PM[k]
            if abs(SATX[I]) < ZERO:
                ALG[I] = 1.0E10
            else:
                ALG[I] = (PM[k] - abs(CM[I])) / abs(SATX[I])
        
        # Find minimum load factor - vectorized
        SALG = 1.0E10
        NPH = 0
        TEST_np = CM_np * SATX_np
        ALG_np = ALG._data[1:M2+1]
        
        # Only consider positive tests (moment increasing toward Mp)
        valid_mask = TEST_np >= 0.0
        if np.any(valid_mask):
            valid_ALG = np.where(valid_mask, ALG_np, 1.0E10)
            min_idx = np.argmin(valid_ALG)
            SALG = valid_ALG[min_idx]
            NPH = min_idx + 1
        
        # Scale forces and displacements - vectorized
        SATX_np = SATX._data[1:M3+1] * SALG
        for I in range(1, M3+1):
            SATX[I] = SATX_np[I-1]
        
        CLG += SALG
        
        # Update cumulative values - vectorized
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
        
        # Determine which member and joint
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
        fp.write(f"%     * PLASTIC HINGE {NCYCL:3d} FORMED IN MEMBER {I:3d} NEAR JOINT {J:3d} WHEN LOAD FACTOR IS {CLG:12.3f}\n%\n")
        
        fp.write("%          CUMULATIVE DEFORMATIONS\n")
        fp.write("%               JOINT    X-DISP       Y-DISP       ROTN\n")
        
        LL = 0
        for I in range(1, JCT+1):
            disp_out = [0.0, 0.0, 0.0]
            for J in range(1, 4):
                if JTYPE[I,J]:
                    LL += 1
                    disp_out[J-1] = CX[LL]
            fp.write(f"%{I:19d}{disp_out[0]:13.5f}{disp_out[1]:13.5f}{disp_out[2]:13.5f}\n")
        
        fp.write("%\n%          CUMULATIVE MOMENTS\n")
        fp.write("%               MEMBER       END MOMENTS            JOINTS     PLASTIC MOM\n")
        
        for I in range(1, NM+1):
            k = 2*I - 1
            fp.write(f"%{I:19d}{CM[k]:14.2f}{CM[k+1]:11.2f}{MCON[I,1]:7d} AND{MCON[I,2]:2d}{PM[I]:14.2f}\n")
        
        fp.write("%\n%          CUMULATIVE TENSION FORCES\n")
        fp.write("%               MEMBER     TENSION\n")
        
        for I in range(1, NM+1):
            fp.write(f"%{I:19d}{CT[I]:15.2f}\n")
        
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
    
    fp.write(f"%\n%     ANALYSIS COMPLETED FOR FRAME NO {JFN:3d}\n\n")
    fp.close()
    print(f"Analysis completed. Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python epframe.py input_file output_file")
        sys.exit(1)
    
    epframe_analysis(sys.argv[1], sys.argv[2])
