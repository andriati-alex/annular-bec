#include "system_solvers.h"
#include "array_memory.h"
#include "array_operations.h"

/* Program to test the accuracy Numerical error of LU decomposition
 * ===============================================================
 *
 * Test System taken from example 7.4
 *
 */

int main() {
    double condu;
    double errback;

    Carray upper = carrDef(2);
    Carray lower = carrDef(2);
    Carray mid = carrDef(3);

    upper[0] = sqrt(3.0);
    upper[1] = sqrt(5.0);
    lower[0] = - sqrt(2.0) * (sqrt(5) + 1E-14);
    lower[1] = - sqrt(3.0) * (1 + 1E-12);
    mid[0] = -sqrt(2.0);
    mid[1] = sqrt(3.0) * 1E-14;
    mid[2] = sqrt(5.0) * 1E-12;

    condu = cond(3, upper, lower, mid);
    errback = errBack(3, upper, lower, mid);

    printf("cond = %.3g, errback = %.3g", condu, errback);

    printf("\n");
    return 0;
}
