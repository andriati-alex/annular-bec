#include "../include/rmatrix_operations.h"



/*          ***********************************************          */
/*                   SETUP VALUES IN MATRIX ENTRIES                  */
/*          ***********************************************          */



void rmatFill(int m, int n, double z, Rmatrix M)
{   /* Fill all matrix with a single constant value */
    unsigned int i, j;
    for (i = 0; i < m; i++) { for (j = 0; j < n; j++) { M[i][j] = z; } }
}

void rmatFillDK(int n, int k, Rarray z, Rmatrix M)
{   /* Fill a diagonal starting from k-column or row */
    unsigned int i;
    if (k >= 0) { for (i = 0; i < n - k; i++) M[i][i+k] = z[i]; }
    else        { for (i = 0; i < n + k; i++) M[i-k][i] = z[i]; }
}

void rmatFillTri(int n, Rarray upper, Rarray mid, Rarray lower, Rmatrix M)
{   /* Fill trigiaonal matrix */
    rmatFill(n, n, 0, M);
    rmatFillDK(n, -1, lower, M);
    rmatFillDK(n,  1, upper, M);
    rmatFillDK(n,  0, mid, M);
}

RCCSmat emptyRCCS(int n, int max_nonzeros)
{
    RCCSmat M = (struct RCCS * restrict) malloc(sizeof(struct RCCS));
    M->m = max_nonzeros;
    M->vec = rarrDef(max_nonzeros * n);
    M->col = (int *) malloc(max_nonzeros * n * sizeof(int));

    return M;
}

void setValueRCCS(int n, int i, int j, int col, double z, RCCSmat M)
{   /* j-th non-zero term in the row i and col its true column position */
    M->vec[i + n * j] = z;
    M->col[i + n * j] = col;
}

void addValueRCCS(int n, int i, int j, double z, RCCSmat M)
{ M->vec[i + n * j] += z; }

RCCSmat triToRCCS(int n, Rarray upper, Rarray lower, Rarray mid)
{   /* Given three diagonals construct CCS matrix format */
    unsigned int j;

    RCCSmat M = (struct RCCS * restrict)malloc(sizeof(struct RCCS));

    M->m = 3;
    M->vec = rarrDef(3 * n);
    M->col = (int *) malloc(3 * n * sizeof(int));

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
    for (j = 2 * n + 1; j < 3 * n - 1; j++) { M->col[j] = j + 1 - 2 * n; }

    return M;
}

RCCSmat CyclicToRCCS(int n, Rarray upper, Rarray lower, Rarray mid)
{   /* Given three diagonals and cyclic terms construct CCS format */
    unsigned int j;
    
    RCCSmat M = (struct RCCS * restrict)malloc(sizeof(struct RCCS));

    M->m = 3;
    M->vec = rarrDef(3 * n);
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



/*          **********************************************          */
/*          MATRIX-VECTOR AND MATRIX-MATRIX MULTIPLICATION          */
/*          **********************************************          */



void rmatvec(int m, int n, Rmatrix M, Rarray v, Rarray ans) {
    unsigned int i, j;
    for (i = 0; i < m; i++) {
        ans[i] = M[i][0] * v[0];
        for (j = 1; j < n; j++) { ans[i] += M[i][j] * v[j]; }
    }
}

void rmatmat(int m, int n, int l, Rmatrix M, Rmatrix A, Rmatrix ans)
{
    unsigned int i, j, k;
    for (i = 0; i < m; i++) { 
        for (j = 0; j < l; j++) {
            for (k = 0; k < n; k++) ans[i][j] += M[i][k] * A[k][j];
        }
    }
}

void RCCSvec(int n, Rarray vals, int * restrict cols, int m,
             Rarray vec, Rarray ans)
{
    unsigned int i, l;
    #pragma omp parallel for private(l, i)
    for (i = 0; i < n; i++) {
        ans[i] = vals[i] * vec[cols[i]];
        for (l = 1; l < m; l++) {
            ans[i] += vals[i + l * n] * vec[cols[i + l * n]];
        }
    }
}
