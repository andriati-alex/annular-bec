#include "system_solvers.h"
#include "array_memory.h"
#include "array_operations.h"

/* Program to test the solvers with large systems
 * **********************************************/

int main() {
    
    int l;
    double complex det;

    Carray upper = carrDef(3);
    Carray lower = carrDef(3);
    Carray mid = carrDef(3);
    Carray RHS = carrDef(3);
    Carray ans = carrDef(3);

    Cmatrix A = cmatDef(3, 3);
    Cmatrix M = cmatDef(3, 3);

    CCSmat Accs;

    /**************** Test 1 ****************/
    // Independent of what we set in cyclic position
    // that is the last value of upper(top term) and lower(bottom term)

    printf("\n\n\t================= Test 1 ================= \n\n");

    upper[0] = 1 + 1 * I;
    upper[1] = 1;
    lower[0] = 1 - 1 * I;
    lower[1] = 1;
    mid[0] = 0;
    mid[1] = 3;
    mid[2] = 2;

    RHS[0] =  2;
    RHS[1] =  4 + 2 * I;
    RHS[2] = -4 - 2 * I;

    carrFill(3, 0, ans);
    
    cmatFillTri(3, upper, mid, lower, A);
    Accs = triToCCS(3, upper, lower, mid);

    // Pre-Conditioner
    invTri(upper, lower, mid, 3, M);

    /**************** Conjugate-Gradient ****************/
    printf("\n\n\t *** Conjugate Gradient ***");

    l = CCG(3, A, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient ***");

    l = preCCG(3, A, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(use M) ***");

    l = MpreCCG(3, A, RHS, ans, 1E-8, upper, lower, mid);
    
    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Conjugate Gradient(CCS) ***");

    l = CCSCCG(3, Accs, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);

    /**************** Pre-Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) ***");

    l = preCCSCCG(3, Accs, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(CCS) (use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) (use M) ***");

    l = MpreCCSCCG(3, Accs, RHS, ans, 1E-8, upper, lower, mid);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Tridiagonal LU ****************/
    printf("\n\n\t *** Tridiagonal LU ***");

    triDiag(3, upper, lower, mid, RHS, ans);

    det = checkLU(3, upper, lower, mid);
    if (cabs(det) == 0) { printf("\nnull determinant, method failed.\n"); }

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Tridiagonal LU-Cyclic ****************/
    printf("\n\n\t *** Tridiagonal LU-Cyclic(first main diag. = 0)***");

    upper[2] = 0; // simple tridiagonal
    lower[2] = 0;
    triCyclicLU(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);

    carrFill(3, 0, ans);

    /**************** Tridiagonal-Cyclic Sherman-Morrison ****************/
    printf("\n\n\t *** Tridiagonal-Cyclic Sherman-Morrison ***");

    upper[2] = 0; // simple tridiagonal
    lower[2] = 0;

    triCyclicSM(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);

    carrFill(3, 0, ans);

    /**************** Test 2 ****************/
    // Independent of what we set in cyclic position
    // that is the last value of upper(top term) and lower(bottom term)
    
    printf("\n\n\n\t================= Test 2 ================= \n\n");

    upper[0] = 1 + 1 * I;
    upper[1] = 2 - 1 * I;
    lower[0] = 1 - 1 * I;
    lower[1] = 2 + 1 * I;
    mid[0] = 1;
    mid[1] = 2;
    mid[2] = 2;

    RHS[0] = 1;
    RHS[1] = 0.5;
    RHS[2] = 2 * I;

    carrFill(3, 0, ans);

    cmatFillTri(3, upper, mid, lower, A);
    Accs = triToCCS(3, upper, lower, mid);

    // Pre-Conditioner
    invTri(upper, lower, mid, 3, M);


    /**************** Conjugate-Gradient ****************/
    printf("\n\n\t *** Conjugate Gradient ***");

    l = CCG(3, A, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient ***");

    l = preCCG(3, A, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(use M) ***");

    l = MpreCCG(3, A, RHS, ans, 1E-8, upper, lower, mid);
    
    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Conjugate Gradient(CCS) ***");

    l = CCSCCG(3, Accs, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);

    /**************** Pre-Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) ***");

    l = preCCSCCG(3, Accs, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(CCS) (use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) (use M) ***");

    l = MpreCCSCCG(3, Accs, RHS, ans, 1E-8, upper, lower, mid);

    carrPrint(3, ans);

    carrFill(3, 0, ans);

    /**************** Tridiagonal LU ****************/
    printf("\n\n\t *** Tridiagonal LU ***");

    triDiag(3, upper, lower, mid, RHS, ans);

    det = checkLU(3, upper, lower, mid);
    if (cabs(det) == 0) { printf("\nnull determinant, method failed.\n"); }

    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Tridiagonal LU-Cyclic ****************/
    printf("\n\n\t *** Tridiagonal LU-Cyclic ***");

    upper[2] = 0; // simple tridiagonal
    lower[2] = 0;
    triCyclicLU(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);
    
    carrFill(3, 0, ans);

    /**************** Tridiagonal-Cyclic Sherman-Morrison ****************/
    printf("\n\n\t *** Tridiagonal-Cyclic Sherman-Morrison ***");
    
    upper[2] = 0; // Simple Tridiagonal
    lower[2] = 0;
    
    triCyclicSM(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Test 3 ****************/
    // For now we get cyclic matrices so if we try to apply
    // LU decomposition it's going to solve without the cyclic terms
    
    printf("\n\n\n\t================= Test 3 ================= \n\n");

    upper[0] = 1 + 1 * I;
    upper[1] = 2 - 1 * I;
    lower[0] = 1 - 1 * I;
    lower[1] = 2 + 1 * I;
    mid[0] = 0;
    mid[1] = 3;
    mid[2] = 2;

    RHS[0] =  2;
    RHS[1] =  4 + 2 * I;
    RHS[2] = -4 - 2 * I;

    // Cyclic terms
    upper[2] =  1 * I;
    lower[2] = -1 * I;

    carrFill(3, 0, ans);

    cmatFillTri(3, upper, mid, lower, A);
    A[0][2] = upper[2];
    A[2][0] = lower[2];

    Accs = CyclicToCCS(3, upper, lower, mid);

    // Pre-Conditioner
    invTri(upper, lower, mid, 3, M);

    /**************** Conjugate-Gradient ****************/
    printf("\n\n\t *** Conjugate Gradient ***");

    l = CCG(3, A, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient ***");

    l = preCCG(3, A, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(use M) ***");

    l = MpreCCG(3, A, RHS, ans, 1E-8, upper, lower, mid);
    
    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Conjugate Gradient(CCS) ***");

    l = CCSCCG(3, Accs, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);

    /**************** Pre-Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) ***");

    l = preCCSCCG(3, Accs, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(CCS) (use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) (use M) ***");

    l = MpreCCSCCG(3, Accs, RHS, ans, 1E-8, upper, lower, mid);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Tridiagonal LU ****************/
    printf("\n\n\t *** Tridiagonal LU ***");

    triDiag(3, upper, lower, mid, RHS, ans);

    det = checkLU(3, upper, lower, mid);
    if (cabs(det) == 0) { printf("\nnull determinant, method failed.\n"); }

    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Tridiagonal LU-Cyclic ****************/
    printf("\n\n\t *** Tridiagonal LU-Cyclic(first main diag. = 0)***");

    triCyclicLU(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);
    
    carrFill(3, 0, ans);

    /**************** Tridiagonal-Cyclic Sherman-Morrison ****************/
    printf("\n\n\t *** Tridiagonal-Cyclic Sherman-Morrison ***");
    
    triCyclicSM(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Test 4 ****************/
    // For now we get cyclic matrices so if we try to apply
    // LU decomposition it's going to solve without the cyclic terms
    
    printf("\n\n\n\t================= Test 4 ================= \n\n");

    upper[0] = 1 + 1 * I;
    upper[1] = 1;
    lower[0] = 1 - 1 * I;
    lower[1] = 1;
    mid[0] = -1;
    mid[1] =  3;
    mid[2] =  2;

    RHS[0] =  2;
    RHS[1] =  4 + 2 * I;
    RHS[2] = -4 - 2 * I;

    // Cyclic terms
    upper[2] =  1 * I;
    lower[2] = -1 * I;

    carrFill(3, 0, ans);

    cmatFillTri(3, upper, mid, lower, A);
    A[0][2] = upper[2];
    A[2][0] = lower[2];

    Accs = CyclicToCCS(3, upper, lower, mid);

    // Pre-Conditioner
    invTri(upper, lower, mid, 3, M);

    /**************** Conjugate-Gradient ****************/
    printf("\n\n\t *** Conjugate Gradient ***");

    l = CCG(3, A, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient ***");

    l = preCCG(3, A, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(use M) ***");

    l = MpreCCG(3, A, RHS, ans, 1E-8, upper, lower, mid);
    
    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Conjugate Gradient(CCS) ***");

    l = CCSCCG(3, Accs, RHS, ans, 1E-8);

    carrPrint(3, ans);

    carrFill(3, 0, ans);

    /**************** Pre-Conjugate-Gradient(CCS) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) ***");

    l = preCCSCCG(3, Accs, RHS, ans, 1E-8, M);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Pre-Conjugate-Gradient(CCS) (use M) ****************/
    printf("\n\n\t *** Pre-Conjugate Gradient(CCS) (use M) ***");

    l = MpreCCSCCG(3, Accs, RHS, ans, 1E-8, upper, lower, mid);

    carrPrint(3, ans);

    carrFill(3, 0, ans);
    
    /**************** Tridiagonal LU ****************/
    printf("\n\n\t *** Tridiagonal LU ***");

    triDiag(3, upper, lower, mid, RHS, ans);

    det = checkLU(3, upper, lower, mid);
    if (cabs(det) == 0) { printf("\nnull determinant, method failed.\n"); }

    carrPrint(3, ans);
    
    carrFill(3, 0, ans);
    
    /**************** Tridiagonal LU-Cyclic ****************/
    printf("\n\n\t *** Tridiagonal LU-Cyclic(first main diag. = 0)***");

    triCyclicLU(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);
    
    carrFill(3, 0, ans);

    /**************** Tridiagonal-Cyclic Sherman-Morrison ****************/
    printf("\n\n\t *** Tridiagonal-Cyclic Sherman-Morrison ***");
    
    triCyclicSM(3, upper, lower, mid, RHS, ans);

    carrPrint(3, ans);

    printf("\n");
    return 0;
}
