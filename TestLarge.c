#include "system_solvers.h"
#include "array_memory.h"
#include "array_operations.h"

#include <time.h>
// For measure time see:
// https://www.geeksforgeeks.org/how-to-measure-time-taken-by-a-program-in-c/

#define PI 3.1415

/* Program to test the solvers for Large systems */
/* ********************************************* */

int main() {
    
    int n = 20000;
    int i, l;

    clock_t start, end;
    double cpu_time_used;

    Carray upper = carrDef(n);
    Carray lower = carrDef(n);
    Carray mid = carrDef(n);
    Carray RHS = carrDef(n);
    Carray ans = carrDef(n);

    carrFill(n, 0.79 + 0.63 * I, upper);
    carrFill(n, 0.79 - 0.63 * I, lower);
    carrFill(n, 0.35, mid);
    carrFill(n, 0.0, ans);
    for (i = 0; i < n; i++) { 
        RHS[i] = 0.5 + cos(10  * i * PI / n) + I * sin(4 * I * PI / n);
    }

    // Compressed Column format
    start = clock();
    CCSmat Accs = CyclicToCCS(n, upper, lower, mid);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime taken for CCS storage: %.4f", cpu_time_used);

    /**************** Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Conjugate Gradient(CCS) ***");

    start = clock();
    l = CCSCCG(n, Accs, RHS, ans, 1E-8);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.4f, iterated %d times", cpu_time_used, l);

    carrPrint(n, ans);

    carrFill(n, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(CCS) (use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) (use M) ***");

    start = clock();
    l = MpreCCSCCG(n, Accs, RHS, ans, 1E-8, upper, lower, mid);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.4f, iterated %d times", cpu_time_used, l);

    carrPrint(n, ans);

    carrFill(n, 0, ans);
    
    /**************** Tridiagonal LU ****************/
    printf("\n\n\t *** Tridiagonal LU ***");

    start = clock();
    triDiag(n, upper, lower, mid, RHS, ans);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.4f", cpu_time_used);

    carrPrint(n, ans);

    carrFill(n, 0, ans);
    
    /**************** Tridiagonal LU-Cyclic ****************/
    printf("\n\n\t *** Tridiagonal LU-Cyclic(first main diag. = 0)***");

    start = clock();
    triCyclicLU(n, upper, lower, mid, RHS, ans);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.4f", cpu_time_used);

    carrPrint(n, ans);

    carrFill(n, 0, ans);

    /**************** Tridiagonal-Cyclic Sherman-Morrison ****************/
    printf("\n\n\t *** Tridiagonal-Cyclic Sherman-Morrison ***");

    start = clock();
    triCyclicSM(n, upper, lower, mid, RHS, ans);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.4f", cpu_time_used);

    carrPrint(n, ans);

    printf("\n");
    return 0;
}
