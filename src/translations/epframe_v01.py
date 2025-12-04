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
    LU decomposition with partial pivoting
    Returns L, U, and permutation vector
    """
    L = Index1Array((n, n))
    U = Index1Array((n, n))
    perm = Index1Array(n, dtype=int)
    
    # Copy A to U
    for i in range(1, n+1):
        perm[i] = i
        for j in range(1, n+1):
            U[i,j] = A[i,j]
    
    for k in range(1, n+1):
        # Find pivot
        max_val = abs(U[k,k])
        pivot = k
        for i in range(k+1, n+1):
            if abs(U[i,k]) > max_val:
                max_val = abs(U[i,k])
                pivot = i
        
        # Swap rows if necessary
        if pivot != k:
            for j in range(1, n+1):
                U[k,j], U[pivot,j] = U[pivot,j], U[k,j]
                if j < k:
                    L[k,j], L[pivot,j] = L[pivot,j], L[k,j]
            perm[k], perm[pivot] = perm[pivot], perm[k]
        
        # Check for singular matrix
        if abs(U[k,k]) < 1e-20:
            return None, None, None
        
        # Elimination
        for i in range(k+1, n+1):
            L[i,k] = U[i,k] / U[k,k]
            for j in range(k, n+1):
                U[i,j] = U[i,j] - L[i,k] * U[k,j]
    
    # Set diagonal of L to 1
    for i in range(1, n+1):
        L[i,i] = 1.0
    
    return L, U, perm

def lu_solve(A, b, n):
    """
    Solve Ax = b using LU decomposition
    Returns solution vector x (1-indexed)
    """
    L, U, perm = lu_decompose(A, n)
    if L is None:
        return None
    
    # Forward substitution: Ly = Pb
    y = Index1Array(n)
    for i in range(1, n+1):
        sum_val = 0.0
        for j in range(1, i):
            sum_val += L[i,j] * y[j]
        y[i] = b[perm[i]] - sum_val
    
    # Back substitution: Ux = y
    x = Index1Array(n)
    for i in range(n, 0, -1):
        sum_val = 0.0
        for j in range(i+1, n+1):
            sum_val += U[i,j] * x[j]
        x[i] = (y[i] - sum_val) / U[i,i]
    
    return x

def read_input_file(filename):
    """Read input file and return problem data"""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    idx = 0
    JFN = int(lines[idx].split()[0])
    idx += 1
    
    parts = lines[idx].split()
    JCT, NM, E = int(parts[0]), int(parts[1]), float(parts[2])
    idx += 1
    
    # Joint data
    CORD = Index1Array((JCT, 2))
    JTYPE = Index1Array((JCT, 3), dtype=int)
    
    for i in range(1, JCT+1):
        parts = lines[idx].split()
        joint_num = int(parts[0])
        CORD[i,1] = float(parts[1])
        CORD[i,2] = float(parts[2])
        JTYPE[i,1] = int(parts[3])
        JTYPE[i,2] = int(parts[4])
        JTYPE[i,3] = int(parts[5])
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
    
    # Input Echo: Joint Data
    fp.write("%\n%     * DATA FOR JOINTS\n")
    fp.write("%          JOINT   X-COORD   Y-COORD   DFX   DFY   DFZ\n%\n")
    for i in range(1, JCT+1):
        fp.write(f"%{i:14d}{CORD[i,1]:12.2f}{CORD[i,2]:10.2f}{JTYPE[i,1]:6d}{JTYPE[i,2]:6d}{JTYPE[i,3]:6d}\n")
    
    # Input Echo: Member Data
    fp.write("%\n%     * DATA FOR MEMBERS\n")
    fp.write("%          MEMBER   JT1     JT2       IXX      AREA        MP\n%\n")
    for i in range(1, NM+1):
        fp.write(f"%{i:14d}{MCON[i,1]:9d}{MCON[i,2]:8d}{SMA[i]:10.2f}{AREA[i]:10.2f}{PM[i]:10.2f}\n")
    
    # Input Echo: Load Data
    fp.write("%\n%     * DATA FOR LOADS\n")
    fp.write("%          JOINT        PX        PY        PZ\n")
    for load in loads:
        JN, FX, FY, FZ = load
        fp.write(f"%{JN:14d}{FX:12.2f}{FY:10.2f}{FZ:10.2f}\n")
    
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
    
    # Write initial state record (cycle 0)
    # Format: NCYCL, member, joint, CLG, then CX values, CM values, CT values
    fp.write("%\n")
    fp.write(f"{0:4d}{0:4d}{0:4d}")
    fp.write(f"{0.0:16.4E}")  # CLG = 0
    # CX values (all zeros) - need 3*JCT slots but only L are active
    for i in range(1, 3*JCT+1):
        fp.write(f"{0.0:16.4E}")
    # CM values (2*NM)
    for i in range(1, M2+1):
        fp.write(f"{CM[i]:16.4E}")
    # CT values (NM)
    for i in range(1, NM+1):
        fp.write(f"{CT[i]:16.4E}")
    fp.write("\n")
    
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
    
    # Main analysis loop
    while True:
        NCYCL += 1
        
        # Form stiffness matrix KSAT
        CSAT = Index1Array(M3)
        KSAT = Index1Array((L, L+1))
        
        for J in range(1, L+1):
            # Multiply stiffness by compatibility
            for I in range(1, M2+1):
                k = ((I+1)//2)*2 - 1
                CSAT[I] = SF[I,1]*K[J,k] + SF[I,2]*K[J,k+1]
            
            for I in range(1, NM+1):
                k = M2 + I
                CSAT[k] = SA[I] * K[J,k]
            
            # Form KSAT = K^T * S * K
            for I in range(1, L+1):
                KSAT[I,J] = 0.0
                for Kk in range(1, M3+1):
                    KSAT[I,J] += K[I,Kk] * CSAT[Kk]
        
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
        max_disp = 0.0
        for I in range(1, L+1):
            if abs(disp[I]) > max_disp:
                max_disp = abs(disp[I])
        
        if max_disp > XLMT:
            fp.write(f"%\n%     *** DEFORMATIONS LARGER THAN {XLMT:.1f} IN CYCLE NO {NCYCL:4d}\n%\n")
            break
        
        # Calculate member forces
        for I in range(1, M3+1):
            CSAT[I] = 0.0
            for J in range(1, L+1):
                CSAT[I] += K[J,I] * disp[J]
        
        SATX = Index1Array(M3)
        for I in range(1, M2+1):
            k = ((I+1)//2)*2 - 1
            SATX[I] = SF[I,1]*CSAT[k] + SF[I,2]*CSAT[k+1]
        
        for I in range(1, NM+1):
            k = M2 + I
            SATX[k] = SA[I] * CSAT[k]
        
        # Calculate load factors to plastic hinge
        ALG = Index1Array(M2)
        for I in range(1, M2+1):
            k = (I+1)//2
            ZERO = 0.001 * PM[k]
            if abs(SATX[I]) < ZERO:
                ALG[I] = 1.0E10
            else:
                ALG[I] = (PM[k] - abs(CM[I])) / abs(SATX[I])
        
        # Find minimum load factor
        SALG = 1.0E10
        NPH = 0
        for I in range(1, M2+1):
            TEST = CM[I] * SATX[I]
            if TEST >= 0.0:
                if ALG[I] < SALG:
                    SALG = ALG[I]
                    NPH = I
        
        # Scale forces and displacements
        for I in range(1, M3+1):
            SATX[I] = SALG * SATX[I]
        
        CLG += SALG
        
        for I in range(1, M2+1):
            CM[I] += SATX[I]
        
        for I in range(1, NM+1):
            k = M2 + I
            CT[I] += SATX[k]
        
        for I in range(1, L+1):
            disp[I] *= SALG
            CX[I] += disp[I]
        
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
        
        # Write compact data record for post-processing
        # Format: NCYCL, member, joint, CLG, then CX values, CM values, CT values
        fp.write("%\n")
        fp.write(f"{NCYCL:4d}{II:4d}{JJ:4d}")
        fp.write(f"{CLG:16.4E}")
        # CX values - write for all 3*JCT slots (zeros for restrained DOFs)
        LL = 0
        for i in range(1, JCT+1):
            for j in range(1, 4):
                if JTYPE[i,j]:
                    LL += 1
                    fp.write(f"{CX[LL]:16.4E}")
                else:
                    fp.write(f"{0.0:16.4E}")
        # CM values (2*NM)
        for i in range(1, M2+1):
            fp.write(f"{CM[i]:16.4E}")
        # CT values (NM)
        for i in range(1, NM+1):
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
    
    fp.write(f"%\n%     ANALYSIS COMPLETED FOR FRAME NO {JFN:3d}\n\n")
    fp.close()
    print(f"Analysis completed. Results written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python epframe.py input_file output_file")
        sys.exit(1)
    
    epframe_analysis(sys.argv[1], sys.argv[2])
