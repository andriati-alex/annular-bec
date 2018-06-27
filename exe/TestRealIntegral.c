#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../include/calculus.h"

// result by python 1.0724738536687413

/* Compilation
 *
 * gcc -o TestRealIntegral exe/TestRealIntegral.c 
 * -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 
 * -lmkl_gnu_thread -lmkl_core -lm -fopenmp -L./lib -I./include -lgp -O3
 *
 */

double fun(double x) { return pow(cos(x / PI), 2) * log(1 + x); }

int main() {
    double * restrict f = (double * ) malloc(201 * sizeof(double));
    double * restrict x = (double * ) malloc(201 * sizeof(double));
    double dx = 0.01;

    for (int i = 0; i < 201; i++) {
        x[i] = 0 + i * dx;
        f[i] = fun(x[i]);
    }

    printf("\n\t%.15f", Rsimps(201, f, dx));

    printf("\n");
    return 0;
}
