/***** Header file ******/
#include "matrix_operations.h"

void cmatFill(unsigned int m, unsigned int n, double complex z, Cmatrix M) {
    unsigned int i, j;
    for (i = 0; i < m; i++) { for (j = 0; j < n; j++) { M[i][j] = z; } }
}

void cmatvec(unsigned int m, unsigned int n, Cmatrix M, Carray v, Carray ans) {
    unsigned int i, j;
    for (i = 0; i < m; i++) {
        ans[i] = M[i][0] * v[0];
        for (j = 1; j < n; j++) { ans[i] += M[i][j] * v[j]; }
    }
}

void cmatmat(unsigned int m, unsigned int n, unsigned int l,
             Cmatrix M, Cmatrix A, Cmatrix ans)
{
    unsigned int i, j, k;
    for (i = 0; i < m; i++) { 
        for (j = 0; j < l; j++) {
            for (k = 0; k < n; k++) { ans[i][j] += M[i][k] * A[k][j]; }
        }
    }
}

void cmatFillDK(unsigned int n, int k, Carray z, Cmatrix M) {
    unsigned int i;
    if (k >= 0) { for (i = 0; i < n - k; i++) { M[i][i+k] = z[i]; } }
    else        { for (i = 0; i < n + k; i++) { M[i-k][i] = z[i]; } }
}

void cmatFillTri(unsigned int n, Carray upper, Carray mid, Carray lower,
                 Cmatrix M)
{
    cmatFill(n, n, 0, M);
    cmatFillDK(n, -1, lower, M);
    cmatFillDK(n,  1, upper, M);
    cmatFillDK(n,  0, mid, M);
}

struct CCS *triToCCS(unsigned int n, Carray upper, Carray lower, Carray mid)
{
    unsigned int j;
    
    CCSmat M = (struct CCS *)malloc(sizeof(struct CCS));

    M->m = 3;
    M->vec = carrDef(3 * n);
    M->col = (int *)malloc(3 * n * sizeof(int));

    M->vec[0]         = mid[0];
    M->vec[n]         = upper[0];
    M->vec[2 * n]     = 0;
    M->vec[3 * n - 1] = 0;

    for (j = 1; j < n; j++)                 { M->vec[j] = lower[j - 1];     }
    for (j = n + 1; j < 2 * n; j++)         { M->vec[j] = mid[j - n];       }
    for (j = 2 * n + 1; j < 3 * n - 1; j++) { M->vec[j] = upper[j - 2 * n]; }
    
    M->col[0]         = 0;
    M->col[n]         = 1;
    M->col[2 * n]     = 0;
    M->col[3 * n - 1] = 0;
    for (j = 1; j < n; j++)                 { M->col[j] = j - 1; }
    for (j = n + 1; j < 2 * n; j++)         { M->col[j] = j - n; }
    for (j = 2 * n + 1; j < 3 * n - 1; j++) { M->col[j] = j + 1 - 2 * n;   }

    return M;
}

struct CCS *CyclicToCCS(unsigned int n, Carray upper, Carray lower, Carray mid)
{
    unsigned int j;
    
    CCSmat M = (struct CCS *)malloc(sizeof(struct CCS));

    M->m = 3;
    M->vec = carrDef(3 * n);
    M->col = (int *)malloc(3 * n * sizeof(int));

    M->vec[0]         = mid[0];
    M->vec[n]         = upper[0];
    M->vec[2 * n]     = upper[n-1]; // Cyclic term
    M->vec[n - 1]     = lower[n-1]; // Cyclic term
    M->vec[2*n - 1]   = lower[n-2];
    M->vec[3 * n - 1] = mid[n-1];

    // In the last and first line we need to take of bottom cyclic term
    for (j = 1; j < n - 1; j++)             { M->vec[j] = lower[j - 1];     }
    for (j = n + 1; j < 2 * n - 1; j++)     { M->vec[j] = mid[j - n];       }
    for (j = 2 * n + 1; j < 3 * n - 1; j++) { M->vec[j] = upper[j - 2 * n]; }
    
    
    M->col[0]         = 0;
    M->col[n]         = 1;
    M->col[2 * n]     = n - 1; // Cyclic term
    M->col[n - 1]     = 0;     // Cyclic Term
    M->col[2 * n - 1] = n - 2;
    M->col[3 * n - 1] = n - 1;
    for (j = 1; j < n - 1; j++)             { M->col[j] = j - 1;         }
    for (j = n + 1; j < 2 * n - 1; j++)     { M->col[j] = j - n;         }
    for (j = 2 * n + 1; j < 3 * n - 1; j++) { M->col[j] = j + 1 - 2 * n; }

    return M;
}

void CCSvec(unsigned int n, Carray vals, int * cols, int m, 
            Carray vec, Carray ans)
{
    unsigned int i, l, offset;
    for (l = 0; l < m; l++) {
        offset = n * l;
        for (i = offset; i < n + offset; i++) {
            ans[i - offset] += vals[i] * vec[cols[i]];
        }
    }
}

void cmatvecTri(unsigned int n, Cmatrix M, Carray v, Carray ans) {
    unsigned int i;
    ans[0] = M[0][0] * v[0] + M[0][1] * v[1];
    ans[n-1] = M[n-1][n-2] * v[n-2] + M[n-1][n-1] * v[n-1];
    for (i = 1; i < n - 1; i++) {
        ans[i] = M[i][i] * v[i] + M[i][i-1] * v[i-1] + M[i][i+1] * v[i+1];
    }
}

void cmatvecDiag(unsigned int n, int *U, int *L, Cmatrix M, 
                 Carray v, Carray ans) {
    int i, j, k;
    for (j = 0; j < n; j++) {
        if (U[j] < 0) { break; }
        k = U[j];
        for (i = 0; i < n - k; i++) { ans[i] += M[i][i+k] * v[i+k]; }
    }
    for (j = 0; j < n; j++) {
        if (L[j] < 0) { break; }
        k = L[j];
        for (i = k; i < n; i++) { ans[i] += M[i][i-k] * v[i-k]; }
    }
}

double complex detTrik(Carray upper, Carray lower, Carray mid, unsigned int k)
{
    double complex det0 = 1;
    double complex det1 = mid[0];
    double complex detK;
    if (k == -1) { return 0; }
    if (k ==  0) { return 1; }
    return mid[k-1] * detTrik(upper, lower, mid, k - 1) - \
           lower[k-2] * upper[k-2] * detTrik(upper, lower, mid, k - 2);
}

Carray detTri(Carray upper, Carray lower, Carray mid, unsigned int n)
{
    int i;
    Carray theta = carrDef(n + 1);

    theta[0] = 1;
    theta[1] = mid[0];
    for (i = 1; i < n; i++) {
        theta[i+1] = mid[i] * theta[i] - lower[i-1] * upper[i-1] * theta[i-1];
    }
    return theta;
}

Carray phiTri(Carray upper, Carray lower, Carray mid, unsigned int n)
{
    int i;
    Carray phi = carrDef(n + 1);

    phi[n]   = 1;
    phi[n-1] = mid[n-1];
    for (i = n - 2; i >= 0; i--) {
        phi[i] = mid[i] * phi[i+1] - lower[i] * upper[i] * phi[i+2];
    }
    return phi;
}

void invTri(Carray upper, Carray lower, Carray mid, 
            unsigned int n, Cmatrix ans)
{
    int i, j;

    Carray phi = phiTri(upper, lower, mid, n);
    Carray theta = detTri(upper, lower, mid, n);

    double complex thetaN = theta[n];
    double complex prodU, prodL;

    for (i = 0; i < n; i++) {
        ans[i][i] = theta[i] * phi[i+1] / thetaN;
        prodU = 1;
        prodL = 1;
        for (j = i - 1; j >= 0; j--) {
            prodU *= (-1) * upper[j];
            prodL *= (-1) * lower[j];
            ans[j][i] = prodU * theta[j] * phi[i+1] / thetaN;
            ans[i][j] = prodL * theta[j] * phi[i+1] / thetaN;
        }
    }

    free(phi);
    free(theta);
}
