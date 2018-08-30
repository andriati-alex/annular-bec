/*
 * compile with intel mkl installation:
 *
gcc -o test exe/lapack_inversion_test.c -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lm -lgp -L./lib -I./include
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>
#include "../include/matrix_operations.h"

#define PI 3.141592653589793 // to define matrix entries

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

int main(int argc, char * argv[])
{
    int j, // counter
        i, // counter
        N, // The size of the Matrix
        k;

    double arg;

    sscanf(argv[1], "%d", &N);

    double complex x = 0;

    Cmatrix A = cmatDef(N, N);
    Cmatrix A_inv = cmatDef(N, N);

    for (i = 0; i < N; i++)
    {   // Row major storage
        A[i][i] = 1.5 + sin( 3 * PI * i / N );
        for (j = 0; j < i; j++)
        {
            arg = 10 * PI * ((double) i * j) / (N * N);
            A[i][j] = 2 * sin(arg) + 2 + I * 3 * cos(arg);
            A[j][i] = conj(A[i][j]); // Does not matter the upper tirangular
        }
    }

    i = HermitianInv(N, A, A_inv);

    printf("\n\nLapacke returned : %d\n", i);

    printf("\n\nIf there was any problem print identity: \n");

    for (i = 0; i < N; i++)
    {
        printf("\n\t|");
        for (j = 0; j < N; j++)
        {
            x = 0;
            for (k = 0; k < N; k++) x += A[i][k] * A_inv[k][j];
            printf(" (%6.3lf,%6.3lf) |", creal(x), cimag(x));
        }
    }

    cmatFree(N, A); //if one try to free does not work with N >= 4 !!!!
    cmatFree(N, A_inv);

    printf("\n\n");
    return 0;
}
