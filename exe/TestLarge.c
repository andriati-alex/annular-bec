/* Program to test the solvers for Large systems */
/* ********************************************* */

#include "../include/tridiagonal_solver.h"
#include "../include/iterative_solver.h"

#include <time.h>
// For measure time (Serial Code) see:
// https://www.geeksforgeeks.org/how-to-measure-time-taken-by-a-program-in-c/
// To measure time with parallelized code use omp_get_wtime();

#define PI 3.141592653589793
#define NLOOP 10

// 4096    = 2 ^ 12
// 1048576 = 2 ^ 20


// COMPILE
// Run makefile then
// gcc -o TestLarge TestLarge.c -lm -fopenmp -L./lib -I./include -lgp -O3

int main() {

    /* DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER */
    omp_set_num_threads(omp_get_max_threads() / 2);
    /******************************************************/

    int n = 1048576;
    int i, l;

    double start, end;
    double cpu_time_used;
    double sentinel;
    double complex sentinelC;

    /* To test with complex numbers */
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

    /* To test some real numbers routines */
    
    Rarray rupper = rarrDef(n - 1);
    Rarray rlower = rarrDef(n - 1);
    Rarray rmid = rarrDef(n);
    Rarray rRHS = rarrDef(n);
    Rarray rans = rarrDef(n);

    rarrFill(n - 1, 0.79, rupper);
    rarrFill(n - 1, 0.79, rlower);
    rarrFill(n, 0.35, rmid);
    for (i = 0; i < n; i++) { 
        rRHS[i] = 1.5 + cos(10  * i * PI / n);
    }



    // Compressed Column format
    start  = omp_get_wtime();
    CCSmat Accs = CyclicToCCS(n, upper, lower, mid);
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start);
    printf("\n\n\tTime taken for CCS storage: %.9f", cpu_time_used);



    /**************** Pre-Conjugate-Gradient(CCS) (use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) (use M) ***");

    start = omp_get_wtime();
    for (i = 0; i < NLOOP; i++) {
        carrFill(n, 0.0, ans);
        l = MpreCCSCCG(n, Accs, RHS, ans, 1E-8, upper, lower, mid);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start);
    printf("\n\n\tTime: %.9f, mean of %d runs", cpu_time_used / NLOOP, 100);
    
    carrPrint(n, ans);



    /**************** Tridiagonal LU-Cyclic ****************/
    printf("\n\n\t *** Tridiagonal LU-Cyclic ***");

    start = omp_get_wtime();
    for (i = 0; i < NLOOP; i++) {
        triCyclicLU(n, upper, lower, mid, RHS, ans);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start);
    printf("\n\n\tTime: %.9f", cpu_time_used / NLOOP);
    
    carrPrint(n, ans);



    /**************** Tridiagonal Cyclic Sherman-Morrison ****************/
    printf("\n\n\t *** Tridiagonal Cyclic Sherman-Morrison ***");

    start = omp_get_wtime();
    for (i = 0; i < NLOOP; i++) {
        triCyclicSM(n, upper, lower, mid, RHS, ans);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start);
    printf("\n\n\tTime: %.9f", cpu_time_used / NLOOP);
    
    carrPrint(n, ans);



    /**************** CCS Matrix-Vector multiply ****************/
    printf("\n\n\t *** CCS Matrix-Vector multiply ***");

    start = omp_get_wtime();
    for (i = 0; i < 5 * NLOOP; i++) {
        CCSvec(n, Accs->vec, Accs->col, Accs->m, mid, ans);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start) / (5 * NLOOP);
    printf("\n\n\tTime: %.9f", cpu_time_used);



    /**************** Vector Add Function ****************/
    printf("\n\n\t *** Vector Add Function ***");

    start = omp_get_wtime();
    for (i = 0; i < 5 * NLOOP; i++) {
        carrAdd(n, upper, lower, ans);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start) / (5 * NLOOP);
    printf("\n\n\tTime: %.9f", cpu_time_used);



    /**************** Vector Update Function ****************/
    printf("\n\n\t *** Vector Update Function ***");

    start = omp_get_wtime();
    for (i = 0; i < 5 * NLOOP; i++) {
        carrUpdate(n, upper, 0.3 - 0.4 * I, lower, ans);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start) / (5 * NLOOP);
    printf("\n\n\tTime: %.9f", cpu_time_used);



    /**************** Vector Exp Function ****************/
    printf("\n\n\t *** Vector Exp Function ***");

    start = omp_get_wtime();
    for (i = 0; i < 5 * NLOOP; i++) {
        carrExp(n, 0.75 - 0.43 * I, upper, ans);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start) / (5 * NLOOP);
    printf("\n\n\tTime: %.9f", cpu_time_used);



    /**************** Element-wise Absolute value Function ****************/
    printf("\n\n\t ***  Element-wise abs Function ***");

    start = omp_get_wtime();
    for (i = 0; i < 5 * NLOOP; i++) {
        carrAbs(n, upper, ansR);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start) / (5 * NLOOP);
    printf("\n\n\tTime: %.9f", cpu_time_used);



    /**************** Vector Modulus Function ****************/
    printf("\n\n\t ***  Vector Modulus Function ***");

    start = omp_get_wtime();
    for (i = 0; i < 5 * NLOOP; i++) {
        sentinel = carrMod(n, upper);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start) / (5 * NLOOP);
    printf("\n\n\tTime: %.9f", cpu_time_used);



    /**************** Vector scalar product ****************/
    printf("\n\n\t ***  Vector scalar product ***");

    start = omp_get_wtime();
    for (i = 0; i < 5 * NLOOP; i++) {
        sentinelC = carrDot(n, upper, lower);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start) / (5 * NLOOP);
    printf("\n\n\tTime: %.9f", cpu_time_used);





    /* REAL SYSTEMS */

    printf("\n\n\n\t REAL SYSTEMS ");
    
    // Compressed Column format
    start  = omp_get_wtime();
    RCCSmat rAccs = CyclicToRCCS(n, rupper, rlower, rmid);
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start);
    printf("\n\n\tTime taken for real CCS storage: %.9f", cpu_time_used);
    


    /**************** Tridiagonal Real System ****************/
    printf("\n\n\t *** Tridiagonal Real System ***");

    start = omp_get_wtime();
    for (i = 0; i < NLOOP; i++) {
        rtriDiag(n, rupper, rlower, rmid, rRHS, rans);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start);
    printf("\n\n\tTime: %.9f", cpu_time_used / NLOOP);
    
    rarrPrint(n, rans);
    


    /**************** Pre-Conjugate-Gradient(real CCS) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(real CCS) ***");

    start = omp_get_wtime();
    for (i = 0; i < NLOOP; i++) {
        rarrFill(n, 0.0, rans);
        l = RCG(n, rAccs, rRHS, rans, 1E-8, rupper, rlower, rmid);
    }
    end = omp_get_wtime();
    cpu_time_used = (double) (end - start);
    printf("\n\n\tTime: %.9f, mean of %d runs", cpu_time_used / NLOOP, NLOOP);
    
    rarrPrint(n, rans);



    // Release memory

    free(upper); free(lower); free(mid); free(RHS); free(ans);
    free(rupper); free(rlower); free(rmid); free(rRHS); free(rans);

    CCSFree(Accs);
    RCCSFree(rAccs);

    printf("\n");
    return 0;
}
