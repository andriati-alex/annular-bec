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

int main() {

    /* DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER */
    mkl_set_num_threads(4);
    omp_set_num_threads(4);

    int n = 256;
    int i;

    double dx = (2 * PI) / (n - 1);

    double start, time_used; // show time taken

    // allocate memory
    Carray f = carrDef(n);
    Carray dfdx = carrDef(n);
    Carray dfdxdiff = carrDef(n);
    Rarray freal = rarrDef(n);
    Rarray x = rarrDef(n);
    Carray cints = carrDef(n);
    Rarray rints = rarrDef(n);

    // setup input data
    for (i = 0; i < n; i++) {
        x[i] = -PI + i * dx;
        f[i] = cos(x[i]) + I * sin(x[i]);
        dfdx[i] = f[i];
    }

    carrRPart(n, f, freal);


    /*                        *****************                        */
    /*                        COMPUTE INTEGRALS                        */
    /*                        *****************                        */



    start = omp_get_wtime();
    for (i = 0; i < 1000; i++)
        cints[n - 1] = Functional(n, dx, -1, 0, -1, rints, f);
    time_used = (double) (omp_get_wtime() - start) / 1000;
    printf("\n\n\tFunctional took %.6f ms", time_used * 1000);

    for (i = 0; i < n - 1; i++) {
        cints[i] = Csimps(i + 1, f, dx);
        rints[i] = Rsimps(i + 1, freal, dx);
    }

    start = omp_get_wtime();
    for (i = 0; i < 1000; i++) cints[n - 1] = Csimps(n, f, dx);
    time_used = (double) (omp_get_wtime() - start) / 1000;
    printf("\n\n\tIntegral took %.6f ms", time_used * 1000);



    /*                      *******************                      */
    /*                      COMPUTE DERIVATIVES                      */
    /*                      *******************                      */



    start = omp_get_wtime();
    for (i = 0; i < 1000; i++) dxFFT(n, f, dx, dfdx);
    time_used = ((double) (omp_get_wtime() - start)) / 1000;

    printf("\n\n\tDerivative by FFT took %.6f ms", time_used * 1000);

    start = omp_get_wtime();
    for (i = 0; i < 1000; i++) dxCyclic(n, f, dx, dfdxdiff);
    time_used = ((double) (omp_get_wtime() - start)) / 1000;
    
    printf("\n\n\tDerivative by differences took %.6f ms", time_used * 1000);

    for (i = 0; i < n; i ++)
        cints[i] = conj(f[i]) * dfdxdiff[i];

    printf("\n\n\tFFT - DIFF: %.14lf\n", cimag(Csimps(n, cints, dx)) / 2);



    /*                          ***********                          */
    /*                          RECORD DATA                          */
    /*                          ***********                          */


    free(x);
    free(f);
    free(freal);
    free(dfdx);
    free(dfdxdiff);
    free(cints);
    free(rints);

    printf("\n");
    return 0;
}
