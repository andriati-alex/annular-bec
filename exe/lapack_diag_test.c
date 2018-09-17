/*
 * compile with intel mkl installation:
 *
gcc -o test exe/lapack_inversion_test.c -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lm -lgp -L./lib -I./include
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>

int main(int argc, char * argv[])
{

    int i,
        j,
        ldz,
        n; // size of the system
    
    sscanf(argv[1], "%d", &n);
    ldz = n;

    double
        * d = malloc(n * sizeof(double)),
        * e = malloc(n * sizeof(double)),
        * eigvec = malloc(n * n * sizeof(double));

    for (i = 0; i < n; i++)
    {
        d[i] = 0;
        e[i] = 0;
        for (j = 0; j < n; j++) eigvec[i * n + j] = 0;
    }

    e[0] = 1.0;

    i = LAPACKE_dstev (LAPACK_ROW_MAJOR, 'V', n, d, e, eigvec, ldz);

    if (i == 0) printf("\nSuccess!");

    printf("\n\nEigenvalues: ");
    for (i = 0; i < n; i++) printf("%8.5lf  ", d[i]);

    printf("\nEigenvector: ");
    for (i = 0; i < n; i++)
    {
        printf("\n             ");
        for (j = 0; j < n; j++) printf("%8.5lf  ", eigvec[i*n + j]);
    }

    free(e);
    free(d);
    free(eigvec);

    printf("\n\n");
    return 0;
}
