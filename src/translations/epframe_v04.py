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
    
    # Calculate number of DOFs using numpy sum
    L = int(np.sum(JTYPE.as_np()))
    
    # Initialize load vector (zeros by default from Index1Array)
    VL = Index1Array(L)
    
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
    
    # Initialize cumulative variables (zeros by default from Index1Array)
    CM = Index1Array(M2)
    CT = Index1Array(NM)
    CX = Index1Array(L)
    
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
    
    # Build flexibility matrix K (compatibility matrix) - zeros by default
    K = Index1Array((L, M3))
    
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
        
        # Build the full M3 x M3 member stiffness matrix S_full
        # Structure: [SF_block (M2xM2) | 0      ]
        #            [0                | SA_diag]
        S_full = np.zeros((M3, M3))
        
        # Fill in 2x2 blocks for each member's flexural stiffness
        for i in range(1, NM+1):
            row = 2*i - 2  # 0-indexed
            S_full[row, row] = SF[2*i-1, 1]
            S_full[row, row+1] = SF[2*i-1, 2]
            S_full[row+1, row] = SF[2*i, 1]
            S_full[row+1, row+1] = SF[2*i, 2]
        
        # Fill diagonal for axial stiffness
        for i in range(1, NM+1):
            S_full[M2-1+i, M2-1+i] = SA[i]
        
        # Get 0-indexed numpy views
        K_np = K.as_np()  # L x M3
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
        
        # Calculate member deformations: CSAT = K.T @ disp
        CSAT_np = K_np.T @ disp_np
        
        # Calculate member forces: SATX = S_full @ CSAT
        SATX_np = S_full @ CSAT_np
        
        # Convert back to 1-indexed for remaining calculations
        disp = Index1Array(L)
        disp.from_np(disp_np)
        CSAT = Index1Array(M3)
        CSAT.from_np(CSAT_np)
        SATX = Index1Array(M3)
        SATX.from_np(SATX_np)
        
        # Calculate load factors to plastic hinge using numpy
        # Get numpy views
        SATX_m2 = SATX.as_np()[:M2]  # First M2 elements (moments)
        CM_np = CM.as_np()
        PM_np = PM.as_np()
        
        # For each moment location I, the member is k = (I+1)//2 (1-indexed)
        # In 0-indexed: member_idx = I // 2 for I in 0..M2-1
        member_indices = np.arange(M2) // 2  # 0-indexed member for each moment location
        PM_expanded = PM_np[member_indices]  # Plastic moment for each location
        
        # Calculate ZERO threshold and ALG
        ZERO = 0.001 * PM_expanded
        small_satx = np.abs(SATX_m2) < ZERO
        # Use safe division: replace small values with 1 to avoid divide-by-zero
        safe_denom = np.where(small_satx, 1.0, np.abs(SATX_m2))
        ALG_np = np.where(small_satx, 
                         1.0E10, 
                         (PM_expanded - np.abs(CM_np)) / safe_denom)
        
        # Find minimum load factor where TEST >= 0
        TEST = CM_np * SATX_m2
        valid_mask = TEST >= 0.0
        ALG_masked = np.where(valid_mask, ALG_np, 1.0E10)
        NPH = int(np.argmin(ALG_masked)) + 1  # Convert to 1-indexed
        SALG = ALG_masked[NPH - 1]
        
        # Scale forces and displacements, update cumulative values using numpy
        SATX_np *= SALG
        SATX.from_np(SATX_np)
        
        CLG += SALG
        
        # Update CM (first M2 of SATX)
        CM_np += SATX_np[:M2]
        CM.from_np(CM_np)
        
        # Update CT (last NM of SATX)
        CT.as_np()[:] += SATX_np[M2:]
        
        # Update CX
        disp_np *= SALG
        CX.as_np()[:] += disp_np
        
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
        
        # Calculate and output reactions at support joints
        # Reactions are what the support provides (opposite of member forces on joint)
        fp.write("%\n%          REACTIONS AT SUPPORTS\n")
        fp.write("%               JOINT       RX           RY           MZ\n")
        
        for jt in range(1, JCT+1):
            # Check if this joint has any restrained DOFs (is a support)
            is_support = (JTYPE[jt,1] == 0 or JTYPE[jt,2] == 0 or JTYPE[jt,3] == 0)
            if not is_support:
                continue
            
            Rx, Ry, Rz = 0.0, 0.0, 0.0
            
            # Sum contributions from all members connected to this joint
            for mem in range(1, NM+1):
                J1 = MCON[mem,1]
                J2 = MCON[mem,2]
                
                if J1 != jt and J2 != jt:
                    continue
                
                # Get member geometry (direction from J1 to J2)
                dx = CORD[J2,1] - CORD[J1,1]
                dy = CORD[J2,2] - CORD[J1,2]
                L_mem = sqrt(dx*dx + dy*dy)
                Cx = dx / L_mem  # cosine
                Sy = dy / L_mem  # sine
                
                # Get member forces
                Mm1 = CM[2*mem - 1]  # Moment at J1 end
                Mm2 = CM[2*mem]      # Moment at J2 end
                N = CT[mem]          # Axial force (tension positive)
                V = (Mm1 + Mm2) / L_mem  # Shear (from moment equilibrium)
                
                if J1 == jt:
                    # Member forces on joint at end 1
                    Rx += N * Cx - V * Sy
                    Ry += N * Sy + V * Cx
                    Rz += Mm1
                else:  # J2 == jt
                    # Member forces on joint at end 2
                    Rx += -N * Cx + V * Sy
                    Ry += -N * Sy - V * Cx
                    Rz += Mm2
            
            # Reaction is negative of member forces (what support provides)
            fp.write(f"%{jt:19d}{-Rx:13.2f}{-Ry:13.2f}{-Rz:13.2f}\n")
        
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
