#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <math.h>

#define PI 3.141592653589793 // to compute derivatives with FFT

#define LAPACK_ROW_MAJOR 101

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

    MKL_Complex16 * A = (MKL_Complex16 *) malloc(N * N * sizeof(MKL_Complex16));
    MKL_Complex16 * B = (MKL_Complex16 *) malloc(N * N * sizeof(MKL_Complex16));

    printf("%9.2E\n", -12837.33);

    // Row major storage
    for (i = 0; i < N; i++)
    {
        A[i * N + i].real = 1.5 + sin( 3 * PI * ((double) i) / N );
        B[i * N + i].real = 1;
        B[i * N + i].imag = 0;
        A[i * N + i].imag = 0;
        for (j = 0; j < i; j++)
        {
            arg = 10 * PI * ((double) i * j) / (N * N);
            A[i * N + j].real = sin(arg);
            A[i * N + j].imag = cos(arg);
            // A[j * N + i].real = sin(arg);
            // A[j * N + i].imag = - cos(arg) * sin(arg);

            B[i * N + j].real = 0; // Identity
            B[i * N + j].imag = 0; // Identity
            B[j * N + i].real = 0;
            B[j * N + i].imag = 0;
        }
    }
    
    i = LAPACKE_zhesv(LAPACK_ROW_MAJOR, 'L', N, N, A, N, ipiv, B, N);

    printf("\n\tSuccess: %d, time: %.3lf\n", i, cpu_time_used * 1000);

    for (i = 0; i < N; i++)
    {
        printf("\n\t|");
        for (j = 0; j < N; j++)
        {
            printf("%7.3lf %6.3lf|", B[i * N + j].real, B[i * N + j].imag);
        }
    }

    free(ipiv);

    // free(A);
    free(B);

    printf("\n\n");
    return 0;
}
