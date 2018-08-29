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
        N; // The size of the Matrix

    double arg,
           start,
           cpu_time_used;

    sscanf(argv[1], "%d", &N);

    int * ipiv = (int *) malloc(N * sizeof(int));

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

    printf("\n\nInverted matrix: \n");

    for (i = 0; i < N; i++)
    {
        printf("\n\t|");
        for (j = 0; j < N; j++)
        {
            printf(" (%6.3lf,%6.3lf) |", Id[i * N + j].real, Id[i * N + j].imag);
        }
    }

    // free(ipiv);
    free(A);
    // free(Id);

    printf("\n\n");
    return 0;
}
