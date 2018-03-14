#include "system_solvers.h"
#include "array_memory.h"
#include "array_operations.h"

#include <time.h>
// For measure time see:
// https://www.geeksforgeeks.org/how-to-measure-time-taken-by-a-program-in-c/

#define PI 3.1415
#define NLOOP 10

/* Program to test the solvers for Large systems */
/* ********************************************* */

// 4096    = 2 ^ 12
// 1048576 = 2 ^ 20

int main() {
    
    int n = 1048576;
    int i, l;

    clock_t start, end;
    double cpu_time_used;

    double sentinel;
    double complex sentinelC;

    Rarray ansR = rarrDef(n);

    Carray upper = carrDef(n);
    Carray lower = carrDef(n);
    Carray mid = carrDef(n);
    Carray RHS = carrDef(n);
    Carray ans = carrDef(n);

    carrFill(n, 0.79 + 0.63 * I, upper);
    carrFill(n, 0.79 - 0.63 * I, lower);
    carrFill(n, 0.35, mid);
    for (i = 0; i < n; i++) { 
        RHS[i] = 1.5 + cos(10  * i * PI / n) + I * sin(4 * I * PI / n);
    }

    // Compressed Column format
    start = clock();
    CCSmat Accs = CyclicToCCS(n, upper, lower, mid);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime taken for CCS storage: %.9f", cpu_time_used);

    /**************** Pre-Conjugate-Gradient(CCS) (use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) (use M) ***");

    start = clock();
    for (i = 0; i < NLOOP; i++) {
        carrFill(n, 0.0, ans);
        l = MpreCCSCCG(n, Accs, RHS, ans, 1E-8, upper, lower, mid);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f, mean of %d runs", cpu_time_used / NLOOP, 100);
    
    carrPrint(n, ans);

    /**************** Tridiagonal LU-Cyclic ****************/
    printf("\n\n\t *** Tridiagonal LU-Cyclic ***");

    start = clock();
    for (i = 0; i < NLOOP; i++) {
        triCyclicLU(n, upper, lower, mid, RHS, ans);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f", cpu_time_used / NLOOP);

    carrPrint(n, ans);

    carrFill(n, 0, ans);

    /**************** Tridiagonal-Cyclic Sherman-Morrison ****************/
    printf("\n\n\t *** Tridiagonal-Cyclic Sherman-Morrison ***");

    start = clock();
    for (i = 0; i < NLOOP; i++) {
        triCyclicSM(n, upper, lower, mid, RHS, ans);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f", cpu_time_used / NLOOP);

    carrPrint(n, ans);
    
    /**************** Array Add ****************/
    printf("\n\n\t *** Array Add ***");

    start = clock();
    for (i = 0; i < 5 * NLOOP; i++) {
        carrAdd(n, upper, lower, ans);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f", cpu_time_used / (5 * NLOOP));
    
    /**************** Array Exp ****************/
    printf("\n\n\t *** Array Exp ***");

    start = clock();
    for (i = 0; i < 5 * NLOOP; i++) {
        carrExp(n, 0.5 - 0.43 * I, lower, ans);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f", cpu_time_used / (5 * NLOOP));
    
    /**************** Element-wise Absolute2 ****************/
    printf("\n\n\t *** Element-wise Absolute2 ***");

    start = clock();
    for (i = 0; i < 5 * NLOOP; i++) {
        carrAbs2(n, lower, ansR);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f", cpu_time_used / (5 * NLOOP));

    /**************** Vector Modulus ****************/
    printf("\n\n\t *** Vector Modulus ***");

    start = clock();
    for (i = 0; i < 5 * NLOOP; i++) {
        sentinel = carrMod(n, lower);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f", cpu_time_used / (5 * NLOOP));
    
    /**************** Scalar Product ****************/
    printf("\n\n\t *** Scalar Product ***");

    start = clock();
    for (i = 0; i < 5 * NLOOP; i++) {
        sentinelC = carrDot2(n, lower, upper);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n\n\tTime: %.9f", cpu_time_used / (5 * NLOOP));

    printf("\n");
    return 0;
}
