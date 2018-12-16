#include <stdio.h>
#include "../include/calculus.h"

/* PROGRAM TO TEST DERIVATIVES AND INTEGRALS
 * -----------------------------------------
 *
 * f(x) = cos(2 * PI * x / 4) + I * sin(PI * x / 4 - PI / 2);
 *
 * generate the file calculus_out.dat  that is read in the
 * python script test_calculus.py and plot derivatives and
 * integrals.
 *
 * First run make and then compile:
 * gcc -o TestCalculus src/tests/TestCalculus.c 
 * -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 
 * -lmkl_gnu_thread -lmkl_core -lm -fopenmp -L./lib -I./include -lgp -O3
 *
 * Run the program ./TestCalculus
 *
 * and then the script $ python test_calculus.py
 *
 ******************************************************************************/

int main(int argc, char * argv [])
{

    int
        i,
        n;
    
    sscanf(argv[1], "%d", &n);



    double
        dx,
        rans,
        start,
        time_used;



    double complex
        cans;



    // allocate memory
    Carray
        f = carrDef(n),
        dfdx = carrDef(n);



    Rarray
        freal = rarrDef(n),
        x = rarrDef(n);



    dx = 2 * PI / (n - 1);
    // setup input data to integrate
    for (i = 0; i < n; i++)
    {
        x[i] = -PI + i * dx;
        f[i] = cos(x[i] / 2) + I * sin(x[i] / 4);
        dfdx[i] = f[i];
    }

    carrRPart(n, f, freal);



    /*                        *****************                        */
    /*                        COMPUTE INTEGRALS                        */
    /*                        *****************                        */



    start = omp_get_wtime();
    for (i = 0; i < 1000; i++) cans = Csimps(n, f, dx);
    time_used = (double) (omp_get_wtime() - start) / 1000;
    printf("\n\nIntegral took %.6f ms", time_used * 1000);

    rans = Rsimps(n, freal, dx);
    printf("\n\nf(x) = cos(x / 2) + i sin(x / 4) integrated");
    printf(" in (-PI, PI] give:");
    printf("\n\n\tAnalitically : 4 + i 0");
    printf("\n\n\tSimpson rule : %.10lf + i %.10lf",rans,cimag(cans));
    cans = Ctrapezium(n, f, dx);
    rans = Rtrapezium(n, freal, dx);
    printf("\n\n\tTrapezium    : %.10lf + i %.10lf",rans,cimag(cans));



    /*                      *******************                      */
    /*                      COMPUTE DERIVATIVES                      */
    /*                      *******************                      */



    // setup input data to differentiate

    for (i = 0; i < n; i++)
    {
        x[i] = -PI + i * dx;
        f[i] = cos(x[i]) + I * sin(x[i]);
        dfdx[i] = f[i];
    }

    printf("\n\nf(x) = cos(x) + i sin(x) - f'(-PI) = 0 - i ");



    start = omp_get_wtime();
    for (i = 0; i < 1000; i++) dxFFT(n, f, dx, dfdx);
    time_used = ((double) (omp_get_wtime() - start)) / 1000;

    printf("\n\nDerivative by FFT took %.6f ms", time_used * 1000);
    printf("\n\tf'(-PI) : %.10f + i %.10f",creal(dfdx[0]),cimag(dfdx[0]));



    start = omp_get_wtime();
    for (i = 0; i < 1000; i++) dxFD(n, f, dx, dfdx);
    time_used = ((double) (omp_get_wtime() - start)) / 1000;

    printf("\n\nDerivative by differences took %.6f ms", time_used * 1000);
    printf("\n\tf'(-PI) : %.10f + i %.10f",creal(dfdx[0]),cimag(dfdx[0]));



    free(x);
    free(f);
    free(freal);
    free(dfdx);

    printf("\n\n");
    return 0;
}
