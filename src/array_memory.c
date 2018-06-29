#include "../include/array_memory.h"



/* ======================================================================== */
/* ======================== ALLOCATE MEMORY SECTION ======================= */
/* ======================================================================== */



Rarray rarrDef(int n)
{ return (double * restrict) malloc(n * sizeof(double)); }

Carray carrDef(int n)
{ return (double complex * restrict) malloc(n * sizeof(double complex)); }

Rmatrix rmatDef(int m, int n)
{   /* Allocate real matrix of m rows by n columns */
    int i;
    Rmatrix M;
    M = (double ** restrict) malloc(m * sizeof(double *));
    for (i = 0; i < m; i++) { M[i] = rarrDef(n); }
    return M;
}

Cmatrix cmatDef(int m, int n)
{   /* Allocate complex matrix of m rows by n columns */
    int i;
    Cmatrix M;
    M = (double complex ** restrict) malloc(m * sizeof(double complex *));
    for (i = 0; i < m; i++) { M[i] = carrDef(n); }
    return M;
}



/* ======================================================================== */
/* ======================== RELEASE MEMORY SECTION ======================== */
/* ======================================================================== */



void rarrFree(Rarray v) { free(v); }

void carrFree(Carray v) { free(v); }

void rmatFree(int m, Rmatrix M)
{   /* Release allocated Real Matrix of m rows */
    int i;
    for (i = 0; i < m; i++) { free(M[i]); }
    free(M);
}

void cmatFree(int m, Cmatrix M)
{   /* Release allocated Complex Matrix of m rows */
    int i;
    for (i = 0; i < m; i++) { free(M[i]); }
    free(M);
}

void CCSFree(CCSmat M)
{   /* Release CCS matrix */
    free(M->col);
    free(M->vec);
    free(M);
}

void RCCSFree(RCCSmat M)
{   /* Release CCS matrix */
    free(M->col);
    free(M->vec);
    free(M);
}



/* ======================================================================== */
/* ======================== COMPLEX PRINT SECTION ========================= */
/* ======================================================================== */



void cPrint(double complex z)
{
    if (cimag(z) > 0) { printf("%2.3E +%.3Ei", creal(z), cimag(z)); }
    else { 
        if (cimag(z) == 0) { printf("%2.3E +0.000E+00i", creal(z)); }
        else               { printf("%2.3E %.3Ei", creal(z), cimag(z)); }
    }
}

void carrPrint(int n, Carray v)
{
    int i;
    if (n < 10) {
        printf("\n"); // print all numbers for short arrays
        for (i = 0; i < n; i++) { printf("\n"); cPrint(v[i]); }
    }
    else {
        printf("\n"); // print first and last three elements
        for (i = 0; i < 3; i++)   { printf("\n\t"); cPrint(v[i]); }
        printf("\n\t\t .\n\t\t .\n\t\t .");
        for (i = n; i > n-3; i--) { printf("\n\t"); cPrint(v[i]); }
    }
}

void rarrPrint(int n, Rarray v)
{
    int i;
    if (n < 10) {
        printf("\n");
        for (i = 0; i < n; i++) printf("\n\t%.9f", v[i]);
    }
    else {
        printf("\n");
        for (i = 0; i < 3; i++)   printf("\n\t%.9f", v[i]);
        printf("\n\t\t .\n\t\t .\n\t\t .");
        for (i = n; i > n-3; i--) printf("\n\t%.9f", v[i]);
    }
}

void carr_txt(char fname [], int M, Carray v)
{
    int j;

    FILE * data_file = fopen(fname, "w");

    if (data_file == NULL)
    {   // impossible to open file with the given name
        printf("ERROR: impossible to open file %s\n", fname);
        return;
    }

    for (j = 0; j < M; j ++) {
        if (cimag(v[j]) > 0) {
            fprintf(data_file, "(%.15E+%.15Ej) ", creal(v[j]), cimag(v[j]));
        }
        else {
            if (cimag(v[j]) == 0) {
                fprintf(data_file, "(%.15E+%.15Ej) ", creal(v[j]), 0.0);
            }
            else {
                fprintf(data_file, "(%.15E%.15Ej) ", creal(v[j]), cimag(v[j]));
            }
        }
    }

    fclose(data_file);
}

void cmat_txt(char fname [],
              int N, int row_step, int M, int col_step, Cmatrix S)
{
    int i, j;

    FILE * data_file = fopen(fname, "w");

    if (data_file == NULL)
    {   // impossible to open file with the given name
        printf("ERROR: impossible to open file %s\n", fname);
        return;
    }

    for (i = 0; i < N; i += row_step) {
        for (j = 0; j < M; j += col_step) {
            if (cimag(S[i][j]) > 0) {
                fprintf(data_file, "(%.15E+%.15Ej) ",
                creal(S[i][j]), cimag(S[i][j]));
            }
            else {
                if (cimag(S[i][j]) == 0) {
                    fprintf(data_file, "(%.15E+%.15Ej) ",
                    creal(S[i][j]), 0.0);
                }
                else {
                    fprintf(data_file, "(%.15E%.15Ej) ",
                    creal(S[i][j]), cimag(S[i][j]));
                }
            }
        }
        fprintf(data_file, "\n");
    }
    fclose(data_file);
}
