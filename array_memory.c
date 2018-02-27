/***** Header file *****/
#include "array_memory.h"

// Memory allocation of real vector
Rarray rarrDef(unsigned int n) {
    return (double *) malloc(n * sizeof(double));
}

// Memory allocation of complex vector
Carray carrDef(unsigned int n) {
    return (double complex *) malloc(n * sizeof(double complex));
}

// Memory allocation of real matrix
Rmatrix rmatDef(unsigned int m, unsigned int n) {
    unsigned int i;
    Rmatrix M;
    M = (double **) malloc(m * sizeof(double *));
    for (i = 0; i < m; i++) { M[i] = rarrDef(n); }
    return M;
}

// Memory allocation of complex matrix
Cmatrix cmatDef(unsigned int m, unsigned int n) {
    unsigned int i;
    Cmatrix M;
    M = (double complex **) malloc(m * sizeof(double complex *));
    for (i = 0; i < m; i++) { M[i] = carrDef(n); }
    return M;
}




/* ======================================================================== */
/* ======================== RELEASE MEMORY SECTION ======================== */
/* ======================================================================== */



// Release memory of allocated real vector
void rarrFree(Rarray v) { free(v); }

// Release memory of allocated complex vector
void carrFree(Carray v) { free(v); }

// Release memory of allocated real matrix
void rmatFree(unsigned int m, Rmatrix M) {
    unsigned int i;
    for (i = 0; i < m; i++) { free(M[i]); }
    free(M);
}

// Release memory of allocated complex matrix
void cmatFree(unsigned int m, Cmatrix M) {
    unsigned int i;
    for (i = 0; i < m; i++) { free(M[i]); }
    free(M);
}



/* ======================================================================== */
/* ======================== COMPLEX PRINT SECTION ========================= */
/* ======================================================================== */



void cPrint(double complex z) {
    if (cimag(z) > 0) { printf("%.3f +%.3fi", creal(z), cimag(z)); }
    else { 
        if (cimag(z) == 0) { printf("%.3f +0.000i", creal(z)); }
        else               { printf("%.3f %.3fi", creal(z), cimag(z)); }
    }
}

// print vector as a column matrix
void carrPrint(unsigned int n, Carray v) {
    unsigned int i;
    if (n < 10) {
        printf("\n");
        for (i = 0; i < n; i++) {
            printf("\n");
            cPrint(v[i]);
        }
    }
    else {
        printf("\n");
        for (i = 0; i < 3; i++) {
            printf("\n");
            cPrint(v[i]);
        }
        printf("\t . \n");
        printf("\t . \n");
        printf("\t . \n");
        for (i = n; i > n-3; i--) {
            printf("\n");
            cPrint(v[i]);
        }
    }
}

// print matrix (only for low-dimension matrices)
void cmatPrint(unsigned int m, unsigned int n, Cmatrix M) {
    unsigned int i;
    unsigned int j;
    for (i = 0; i < m; i++) {
        printf("\n|");
        for (j = 0; j < n; j++) {
            printf("\t");
            cPrint(M[i][j]);
        }
        printf("|");
    }
}
