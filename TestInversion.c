#include "system_solvers.h"
#include "array_memory.h"
#include "array_operations.h"

/* Program to test the Inversion of tridiagonal system
 * ***************************************************
 *
 * Result of multiplication must be unit
 *
 * ***************************************************/

int main() {

    Carray upper = carrDef(2);
    Carray lower = carrDef(2);
    Carray mid   = carrDef(3);
    Cmatrix ans  = cmatDef(3, 3);
    Cmatrix A    = cmatDef(3, 3);
    Cmatrix Id   = cmatDef(3, 3);

    upper[0] = 1 + 1 * I;
    upper[1] = 1;
    lower[0] = 1 - 1 * I;
    lower[1] = 1;
    mid[0] = -1;
    mid[1] =  3;
    mid[2] =  2;

    A[0][0] = mid[0];
    A[0][1] = upper[0];
    A[0][2] = 0;
    A[1][0] = lower[0];
    A[1][1] = mid[1];
    A[1][2] = upper[1];
    A[2][0] = 0;
    A[2][1] = lower[1];
    A[2][2] = mid[2];

    invTri(upper, lower, mid, 3, ans);
    cmatmat(3, 3, 3, ans, A, Id);

    cmatPrint(3, 3, Id);

    printf("\n");
    return 0;
}
