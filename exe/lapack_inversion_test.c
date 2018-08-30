/*
 * compile with intel mkl installation:
 *
 * gcc -o test exe/lapack_inversion_test.c  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lm -ldl
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>

#define PI 3.141592653589793 // to define matrix entries

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

int main(int argc, char * argv[])
{
    int j, // counter
        i, // counter
        N, // The size of the Matrix
        k;

    double arg,
           start,
           cpu_time_used;

    sscanf(argv[1], "%d", &N);

    int * ipiv = (int *) malloc(N * sizeof(int));

    MKL_Complex16 x; x.real = 0; x.imag = 0;

    MKL_Complex16 * A = malloc(N * N * sizeof(MKL_Complex16));
    MKL_Complex16 * Acopy = malloc(N * N * sizeof(MKL_Complex16));
    MKL_Complex16 * Id = malloc(N * N * sizeof(MKL_Complex16));

    for (i = 0; i < N; i++)
    {   // Row major storage
        A[i * N + i].real = 1.5 + sin( 3 * PI * i / N );
        A[i * N + i].imag = 0;
        Acopy[i * N + i].real = 1.5 + sin( 3 * PI * i / N );
        Acopy[i * N + i].imag = 0;
        Id[i * N + i].real = 1;
        Id[i * N + i].imag = 0;
        for (j = 0; j < i; j++)
        {
            arg = 10 * PI * ((double) i * j) / (N * N);
            A[i * N + j].real = 2 * sin(arg) + 2;
            A[i * N + j].imag = 3 * cos(arg);
            A[j * N + i].real = 0; // Does not matter the upper tirangular
            A[j * N + i].imag = 0; // when call lapack routine with 'L'
            
            Acopy[i * N + j].real = 2 * sin(arg) + 2;
            Acopy[i * N + j].imag = 3 * cos(arg);
            Acopy[j * N + i].real = Acopy[i * N + j].real;
            Acopy[j * N + i].imag = - Acopy[i * N + j].imag;

            Id[i * N + j].real = 0; // Identity
            Id[i * N + j].imag = 0;
            Id[j * N + i].real = 0;
            Id[j * N + i].imag = 0;
        }
    }
    
    i = LAPACKE_zhesv(LAPACK_ROW_MAJOR, 'L', N, N, A, N, ipiv, Id, N);

    printf("\n\nLapacke returned : %d\n", i);

    printf("\n\nIf there was any problem print identity: \n");

    for (i = 0; i < N; i++)
    {
        printf("\n\t|");
        for (j = 0; j < N; j++)
        {
            x.real = 0;
            x.imag = 0;
            for (k = 0; k < N; k++)
            {
                x.real += Id[i*N + k].real * Acopy[k*N + j].real;
                x.real -= Id[i*N + k].imag * Acopy[k*N + j].imag;
                x.imag += Id[i*N + k].real * Acopy[k*N + j].imag;
                x.imag += Id[i*N + k].imag * Acopy[k*N + j].real;
            }
            printf(" (%6.3lf,%6.3lf) |", x.real, x.imag);
        }
    }

    free(A); //if one try to free does not work with N >= 4 !!!!
    free(Id);
    free(ipiv);

    printf("\n\n");
    return 0;
}
