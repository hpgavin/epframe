/*******************************************************************************
 Program epframe.c  -  ELASTIC PLASTIC ANALYSIS OF A PLANE FRAME
 
 Translated from FORTRAN epframe.f by Hacksoo Lee, Univ. of Michigan, 1986
 C version using NRutil for 1-indexed arrays and LU decomposition
 
 Compile: gcc -o epframe epframe.c NRutil.c -lm
 Usage:   ./epframe input_file output_file
 
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "NRutil.h"

#define PI 3.141592653589793

/* Function prototypes */
void lu_decompose(double **A, int n, int *indx, double *d);
void lu_backsubstitute(double **A, int n, int *indx, double *b);
void read_input_file(char *filename, int *JFN, int *JCT, int *NM, double *E,
                     double **CORD, int **JTYPE, int **MCON, 
                     double *SMA, double *AREA, double *PM,
                     int *LN, int *loaded_joints, double **loads);
void epframe_analysis(char *input_file, char *output_file);

/*-----------------------------------------------------------------------------
LU_DECOMPOSE - LU decomposition with partial pivoting
 Input:  A[1..n][1..n] - matrix to decompose
         n - matrix dimension
 Output: A - contains L and U (L's diagonal is implicitly 1)
         indx[1..n] - row permutation vector
         d - +1 or -1 depending on number of row interchanges
------------------------------------------------------------------------------*/
void lu_decompose(double **A, int n, int *indx, double *d)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double *vv;
    
    vv = dvector(1, n);
    *d = 1.0;
    
    /* Find implicit scaling for each row */
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(A[i][j])) > big) big = temp;
        if (big == 0.0) {
            fprintf(stderr, "Singular matrix in LU decomposition\n");
            exit(1);
        }
        vv[i] = 1.0 / big;
    }
    
    /* Crout's method with partial pivoting */
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = A[i][j];
            for (k = 1; k < i; k++)
                sum -= A[i][k] * A[k][j];
            A[i][j] = sum;
        }
        
        big = 0.0;
        imax = j;
        for (i = j; i <= n; i++) {
            sum = A[i][j];
            for (k = 1; k < j; k++)
                sum -= A[i][k] * A[k][j];
            A[i][j] = sum;
            
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        
        if (j != imax) {
            for (k = 1; k <= n; k++) {
                dum = A[imax][k];
                A[imax][k] = A[j][k];
                A[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        
        indx[j] = imax;
        
        if (A[j][j] == 0.0) A[j][j] = 1.0e-20;
        
        if (j != n) {
            dum = 1.0 / A[j][j];
            for (i = j+1; i <= n; i++)
                A[i][j] *= dum;
        }
    }
    
    free_dvector(vv, 1, n);
}

/*-----------------------------------------------------------------------------
LU_BACKSUBSTITUTE - Solve Ax=b using LU decomposition
 Input:  A[1..n][1..n] - LU decomposed matrix
         n - matrix dimension
         indx[1..n] - permutation vector from lu_decompose
         b[1..n] - right hand side
 Output: b - solution vector
------------------------------------------------------------------------------*/
void lu_backsubstitute(double **A, int n, int *indx, double *b)
{
    int i, ii = 0, ip, j;
    double sum;
    
    /* Forward substitution */
    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i-1; j++) sum -= A[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;
    }
    
    /* Back substitution */
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i+1; j <= n; j++)
            sum -= A[i][j] * b[j];
        b[i] = sum / A[i][i];
    }
}

/*-----------------------------------------------------------------------------
READ_INPUT_FILE - Read input data from file
------------------------------------------------------------------------------*/
void read_input_file(char *filename, int *JFN, int *JCT, int *NM, double *E,
                     double **CORD, int **JTYPE, int **MCON, 
                     double *SMA, double *AREA, double *PM,
                     int *LN, int *loaded_joints, double **loads)
{
    FILE *fp;
    int i, j, joint_num, mem_num;
    
    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open input file '%s'\n", filename);
        exit(1);
    }
    
    /* Read frame number */
    fscanf(fp, "%d", JFN);
    
    /* Read general data */
    fscanf(fp, "%d %d %lf", JCT, NM, E);
    
    /* Read joint data */
    for (i = 1; i <= *JCT; i++) {
        fscanf(fp, "%d %lf %lf %d %d %d",
               &joint_num, &CORD[i][1], &CORD[i][2],
               &JTYPE[i][1], &JTYPE[i][2], &JTYPE[i][3]);
    }
    
    /* Read member data */
    for (i = 1; i <= *NM; i++) {
        fscanf(fp, "%d %d %d %lf %lf %lf",
               &mem_num, &MCON[i][1], &MCON[i][2],
               &SMA[i], &AREA[i], &PM[i]);
    }
    
    /* Read load data */
    fscanf(fp, "%d", LN);
    for (i = 1; i <= *LN; i++) {
        fscanf(fp, "%d %lf %lf %lf",
               &loaded_joints[i], &loads[i][1], &loads[i][2], &loads[i][3]);
    }
    
    fclose(fp);
}

/*-----------------------------------------------------------------------------
EPFRAME_ANALYSIS - Main elastic-plastic frame analysis
------------------------------------------------------------------------------*/
void epframe_analysis(char *input_file, char *output_file)
{
    FILE *fp;
    int JFN, JCT, NM, LN;
    double E;
    double **CORD, **K, **KSAT, **SF, **loads;
    int **JTYPE, **MCON;
    double *SMA, *AREA, *PM, *OLEN, *SA, *VL, *CSAT, *SATX, *ALG;
    double *CM, *CT, *CX, *disp;
    int *loaded_joints, *indx;
    int L, M2, M3, i, j, k, m, J, M, JN;
    int NA, NJ, NK, JF, MJ, MF, NN;
    double X, Y, D, S, C;
    int NCYCL, NPH, I, II, JJ, LL, LJ, ITEST;
    double CLG, SALG, ZERO, TEST, XLMT, max_disp, lu_d;
    double disp_out[3];
    
    /* Allocate temporary arrays for reading */
    CORD = dmatrix(1, 300, 1, 2);
    JTYPE = imatrix(1, 300, 1, 3);
    MCON = imatrix(1, 300, 1, 2);
    SMA = dvector(1, 300);
    AREA = dvector(1, 300);
    PM = dvector(1, 300);
    loaded_joints = ivector(1, 100);
    loads = dmatrix(1, 100, 1, 3);
    
    /* Read input file */
    read_input_file(input_file, &JFN, &JCT, &NM, &E, CORD, JTYPE, MCON,
                    SMA, AREA, PM, &LN, loaded_joints, loads);
    
    /* Open output file */
    if ((fp = fopen(output_file, "w")) == NULL) {
        fprintf(stderr, "Error: cannot open output file '%s'\n", output_file);
        exit(1);
    }
    
    fprintf(fp, "%%\n");
    fprintf(fp, "%%     ELASTIC PLASTIC ANALYSIS OF FRAME NO %3d\n", JFN);
    fprintf(fp, "%%     ---------------------------------------\n%%\n");
    
    /* Calculate number of degrees of freedom */
    L = 0;
    for (i = 1; i <= JCT; i++)
        for (j = 1; j <= 3; j++)
            L += JTYPE[i][j];
    
    /* Allocate working arrays */
    OLEN = dvector(1, NM);
    VL = dvector(1, L);
    SF = dmatrix(1, 2*NM, 1, 2);
    SA = dvector(1, NM);
    K = dmatrix(1, L, 1, 3*NM);
    KSAT = dmatrix(1, L, 1, L+1);
    CSAT = dvector(1, 3*NM);
    SATX = dvector(1, 3*NM);
    ALG = dvector(1, 2*NM);
    CM = dvector(1, 2*NM);
    CT = dvector(1, NM);
    CX = dvector(1, L);
    disp = dvector(1, L);
    indx = ivector(1, L);
    
    /* Initialize load vector */
    for (i = 1; i <= L; i++)
        VL[i] = 0.0;
    
    /* Process loads into load vector */
    for (i = 1; i <= LN; i++) {
        JN = loaded_joints[i];
        LL = 0;
        LJ = JN - 1;
        
        if (LJ > 0) {
            for (j = 1; j <= LJ; j++)
                for (k = 1; k <= 3; k++)
                    LL += JTYPE[j][k];
        }
        
        for (k = 1; k <= 3; k++) {
            if (JTYPE[JN][k]) {
                LL++;
                VL[LL] = loads[i][k];
            }
        }
    }
    
    fprintf(fp, "%%     * GENERAL DATA\n");
    fprintf(fp, "%%          NUMBER OF JOINTS        %6d\n", JCT);
    fprintf(fp, "%%          NUMBER OF MEMBERS       %6d\n", NM);
    fprintf(fp, "%%          MOD OF ELASTICITY  %12.1f\n%%\n", E);
    
    /* Calculate member lengths */
    for (i = 1; i <= NM; i++) {
        int J1 = MCON[i][1];
        int J2 = MCON[i][2];
        X = CORD[J1][1] - CORD[J2][1];
        Y = CORD[J1][2] - CORD[J2][2];
        OLEN[i] = sqrt(X*X + Y*Y);
    }
    
    /* Calculate stiffness coefficients */
    for (i = 1; i <= NM; i++) {
        SF[2*i][2] = 4.0 * E * SMA[i] / OLEN[i];
        SF[2*i-1][1] = SF[2*i][2];
        SF[2*i][1] = 0.5 * SF[2*i][2];
        SF[2*i-1][2] = SF[2*i][1];
        SA[i] = E * AREA[i] / OLEN[i];
    }
    
    M2 = 2 * NM;
    M3 = 3 * NM;
    
    /* Initialize cumulative variables */
    for (i = 1; i <= M2; i++)
        CM[i] = 0.0;
    for (i = 1; i <= NM; i++)
        CT[i] = 0.0;
    for (i = 1; i <= L; i++)
        CX[i] = 0.0;
    
    NCYCL = 0;
    CLG = 0.0;
    
    /* Build compatibility/flexibility matrix K */
    for (i = 1; i <= L; i++)
        for (j = 1; j <= M3; j++)
            K[i][j] = 0.0;
    
    NJ = 0;
    NK = 0;
    
    for (J = 1; J <= JCT; J++) {
        for (M = 1; M <= NM; M++) {
            NA = NJ;
            
            /* Check if joint J is connected to member M */
            if (J == MCON[M][1]) {
                JF = MCON[M][2];
                MJ = 2*M - 1;
                MF = MJ + 1;
            } else if (J == MCON[M][2]) {
                JF = MCON[M][1];
                MJ = 2*M;
                MF = MJ - 1;
            } else {
                continue;
            }
            
            X = CORD[JF][1] - CORD[J][1];
            Y = CORD[JF][2] - CORD[J][2];
            D = sqrt(X*X + Y*Y);
            S = Y / D;
            C = X / D;
            NN = 2*NM + M;
            
            if (JTYPE[J][1]) {
                NA = NA + 1;
                K[NA][MJ] = S / D;
                K[NA][MF] = K[NA][MJ];
                K[NA][NN] = -C;
            }
            
            if (JTYPE[J][2]) {
                NA = NA + 1;
                K[NA][MJ] = -C / D;
                K[NA][MF] = K[NA][MJ];
                K[NA][NN] = -S;
            }
            
            if (JTYPE[J][3]) {
                NA = NA + 1;
                K[NA][MJ] = 1.0;
            }
            
            if (NA > NK)
                NK = NA;
        }
        NJ = NK;
    }
    
    /* Main analysis loop */
    while (1) {
        NCYCL++;
        
        /* Form stiffness matrix KSAT = K^T * S * K */
        for (J = 1; J <= L; J++) {
            /* Multiply stiffness by compatibility: CSAT = S * K[J,:] */
            for (I = 1; I <= M2; I++) {
                k = ((I+1)/2)*2 - 1;
                CSAT[I] = SF[I][1]*K[J][k] + SF[I][2]*K[J][k+1];
            }
            
            for (I = 1; I <= NM; I++) {
                k = M2 + I;
                CSAT[k] = SA[I] * K[J][k];
            }
            
            /* Form KSAT[I][J] = K^T * CSAT */
            for (I = 1; I <= L; I++) {
                KSAT[I][J] = 0.0;
                for (k = 1; k <= M3; k++)
                    KSAT[I][J] += K[I][k] * CSAT[k];
            }
        }
        
        /* Add load vector to augmented matrix */
        for (I = 1; I <= L; I++)
            KSAT[I][L+1] = VL[I];
        
        /* Solve system using LU decomposition */
        lu_decompose(KSAT, L, indx, &lu_d);
        
        for (i = 1; i <= L; i++)
            disp[i] = KSAT[i][L+1];
        
        lu_backsubstitute(KSAT, L, indx, disp);
        
        /* Check for excessive deformations */
        XLMT = 1000.0;
        max_disp = 0.0;
        for (I = 1; I <= L; I++) {
            if (fabs(disp[I]) > max_disp)
                max_disp = fabs(disp[I]);
        }
        
        if (max_disp > XLMT) {
            fprintf(fp, "%%\n%%     *** DEFORMATIONS LARGER THAN %.1f IN CYCLE NO %4d\n%%\n",
                    XLMT, NCYCL);
            break;
        }
        
        /* Calculate member forces: CSAT = K^T * disp */
        for (I = 1; I <= M3; I++) {
            CSAT[I] = 0.0;
            for (J = 1; J <= L; J++)
                CSAT[I] += K[J][I] * disp[J];
        }
        
        /* Calculate moments and axial forces */
        for (I = 1; I <= M2; I++) {
            k = ((I+1)/2)*2 - 1;
            SATX[I] = SF[I][1]*CSAT[k] + SF[I][2]*CSAT[k+1];
        }
        
        for (I = 1; I <= NM; I++) {
            k = M2 + I;
            SATX[k] = SA[I] * CSAT[k];
        }
        
        /* Calculate load factors to plastic hinge formation */
        for (I = 1; I <= M2; I++) {
            k = (I+1)/2;
            ZERO = 0.001 * PM[k];
            if (fabs(SATX[I]) < ZERO) {
                ALG[I] = 1.0E10;
            } else {
                ALG[I] = (PM[k] - fabs(CM[I])) / fabs(SATX[I]);
            }
        }
        
        /* Find minimum positive load factor */
        SALG = 1.0E10;
        NPH = 0;
        for (I = 1; I <= M2; I++) {
            TEST = CM[I] * SATX[I];
            if (TEST >= 0.0) {
                if (ALG[I] < SALG) {
                    SALG = ALG[I];
                    NPH = I;
                }
            }
        }
        
        /* Scale forces and displacements by load factor */
        for (I = 1; I <= M3; I++)
            SATX[I] = SALG * SATX[I];
        
        CLG += SALG;
        
        /* Update cumulative values */
        for (I = 1; I <= M2; I++)
            CM[I] += SATX[I];
        
        for (I = 1; I <= NM; I++) {
            k = M2 + I;
            CT[I] += SATX[k];
        }
        
        for (I = 1; I <= L; I++) {
            disp[I] *= SALG;
            CX[I] += disp[I];
        }
        
        /* Determine which member and joint formed hinge */
        I = (NPH + 1) / 2;
        k = (NPH / 2) * 2 - NPH;
        if (k)
            J = MCON[I][1];
        else
            J = MCON[I][2];
        
        II = I;
        JJ = J;
        
        /* Print results for this cycle */
        fprintf(fp, "%%\n%%\n%%\n");
        fprintf(fp, "%%     * PLASTIC HINGE %3d FORMED IN MEMBER %3d NEAR JOINT %3d WHEN LOAD FACTOR IS %12.3f\n%%\n",
                NCYCL, I, J, CLG);
        
        fprintf(fp, "%%          CUMULATIVE DEFORMATIONS\n");
        fprintf(fp, "%%               JOINT    X-DISP       Y-DISP       ROTN\n");
        
        LL = 0;
        for (I = 1; I <= JCT; I++) {
            for (k = 0; k < 3; k++)
                disp_out[k] = 0.0;
            
            for (J = 1; J <= 3; J++) {
                if (JTYPE[I][J]) {
                    LL++;
                    disp_out[J-1] = CX[LL];
                }
            }
            fprintf(fp, "%%%19d%13.5f%13.5f%13.5f\n", I, disp_out[0], disp_out[1], disp_out[2]);
        }
        
        fprintf(fp, "%%\n%%          CUMULATIVE MOMENTS\n");
        fprintf(fp, "%%               MEMBER       END MOMENTS            JOINTS     PLASTIC MOM\n");
        
        for (I = 1; I <= NM; I++) {
            k = 2*I - 1;
            fprintf(fp, "%%%19d%14.2f%11.2f%7d AND%2d%14.2f\n",
                    I, CM[k], CM[k+1], MCON[I][1], MCON[I][2], PM[I]);
        }
        
        fprintf(fp, "%%\n%%          CUMULATIVE TENSION FORCES\n");
        fprintf(fp, "%%               MEMBER     TENSION\n");
        
        for (I = 1; I <= NM; I++)
            fprintf(fp, "%%%19d%15.2f\n", I, CT[I]);
        
        /* Modify stiffness for plastic hinge */
        ITEST = ((NPH/2)*2) - NPH;
        if (ITEST) {
            SF[NPH+1][2] = 0.75 * SF[NPH+1][2];
            SF[NPH+1][1] = 0.0;
            SF[NPH][1] = 0.0;
            SF[NPH][2] = 0.0;
        } else {
            SF[NPH-1][1] = 0.75 * SF[NPH-1][1];
            SF[NPH-1][2] = 0.0;
            SF[NPH][1] = 0.0;
            SF[NPH][2] = 0.0;
        }
    }
    
    fprintf(fp, "%%\n%%     ANALYSIS COMPLETED FOR FRAME NO %3d\n\n", JFN);
    fclose(fp);
    
    /* Free memory */
    free_dmatrix(CORD, 1, 300, 1, 2);
    free_imatrix(JTYPE, 1, 300, 1, 3);
    free_imatrix(MCON, 1, 300, 1, 2);
    free_dvector(SMA, 1, 300);
    free_dvector(AREA, 1, 300);
    free_dvector(PM, 1, 300);
    free_ivector(loaded_joints, 1, 100);
    free_dmatrix(loads, 1, 100, 1, 3);
    free_dvector(OLEN, 1, NM);
    free_dvector(VL, 1, L);
    free_dmatrix(SF, 1, 2*NM, 1, 2);
    free_dvector(SA, 1, NM);
    free_dmatrix(K, 1, L, 1, 3*NM);
    free_dmatrix(KSAT, 1, L, 1, L+1);
    free_dvector(CSAT, 1, 3*NM);
    free_dvector(SATX, 1, 3*NM);
    free_dvector(ALG, 1, 2*NM);
    free_dvector(CM, 1, 2*NM);
    free_dvector(CT, 1, NM);
    free_dvector(CX, 1, L);
    free_dvector(disp, 1, L);
    free_ivector(indx, 1, L);
    
    printf("Analysis completed. Results written to %s\n", output_file);
}

/*-----------------------------------------------------------------------------
MAIN
------------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s input_file output_file\n", argv[0]);
        exit(1);
    }
    
    epframe_analysis(argv[1], argv[2]);
    
    return 0;
}