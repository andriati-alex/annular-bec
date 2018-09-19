#include "../include/tridiagonal_solver.h"



/*          ***********************************************          */
/*              CHECK SOLVABILITY of LU DECOMPOSITION                */
/*          ***********************************************          */



double cond(int n, Carray upper, Carray lower, Carray mid)
{
    unsigned int i;
    // condition numbers
    double maxCond = 1,
           condK   = 1;

    // L . U decomposition
    double complex u,
                   l;

    u = mid[0];
    condK = 1;
    for (int i = 0;  i < n - 1; i++) {
        l = lower[i] / u;
        u = mid[i+1] - l * upper[i];
        condK = 1 + cabs(upper[i] * l / u) * (2 + condK);
        if (maxCond < condK) { maxCond = condK; }
    }

    return maxCond;
}

double errBack(int n, Carray upper, Carray lower, Carray mid)
{
    unsigned int i;
    // Error backward
    double maxErrBack = 1,
           ErrBack    = 1;

    // L . U decomposition
    double complex u,
                   l;

    u = mid[0];
    for (i = 0;  i < n - 1; i++) {
        l = lower[i] / u;
        u = mid[i+1] - l * upper[i];
        ErrBack = (cabs(l) * cabs(upper[i]) + cabs(u)) / cabs(mid[i+1]);
        if (maxErrBack < ErrBack) { maxErrBack = ErrBack; }
    }

    return maxErrBack * 1E-16;
}

double complex checkLU(int n, Carray upper, Carray lower, Carray mid)
{
    unsigned int i;
    double complex det0 = 1;
    double complex det1 = mid[0];
    double complex detK;
    for (i = 1; i < n; i++) {
        detK = mid[i] * det1 - lower[i-1] * upper[i-1] * det0;
        if (cabs(detK) == 0) { return 0; }
        det0 = det1;
        det1 = detK;
    }
    return detK;
}



/*          ***********************************************          *
 *                        TRIDIAGONAL SOLVERS                        *
 *          ***********************************************          */



void triDiag(int n, Carray upper, Carray lower, Carray mid, Carray RHS,
     Carray ans)
{
    // Intermediate steps - L . U decomposition
    Carray u = carrDef(n);
    Carray l = carrDef(n - 1);
    Carray z = carrDef(n);
    
    // Loop counters variables.
    unsigned int i, k;
    
    if (cabs(mid[0]) == 0 ) {
        // In this case there is a system reduction
        // where we solve for [x1  x3  x4 ... xn]
        // what is equivalent to adjust the two first 
        // equations and starts counters from 1

        double complex RHS1 = RHS[1] - mid[1]   * RHS[0] / upper[0],
                       RHS2 = RHS[2] - lower[1] * RHS[0] / upper[0];

        // u and z factor initizlization with the changed system
        u[1] = lower[0];
        z[1] = RHS1;
        // One iteration need to be performed outer because
        // the change of values in lower and RHS
        l[1] = 0;
        u[2] = mid[2] - l[1] * upper[1];
        z[2] = RHS2 - l[1] * z[1];
        for (i = 2; i < n - 1; i++) {
            k = i + 1;
            l[i] = lower[i] / u[i];
            u[k] = mid[k] - l[i] * upper[i];
            z[k] = RHS[k] - l[i] * z[i];
        }

        ans[n-1] = z[n-1] / u[n-1];
        for (i = 2; i <= n - 1; i++) {
            k = n - i;
            ans[k] = (z[k] - upper[k] * ans[k+1]) / u[k];
        }

        // Obtained order ans[0..n] = [nan  x1  x3  x4 .. xn]
        // Organize ans[0..n] = [x1  x2  x3  .. xn]
        ans[0] = ans[1];
        ans[1] = RHS[0] / upper[0];

        // Free local alilocated memory
        free(u);
        free(l);
        free(z);

        return;
    }
    
    u[0] = mid[0];
    z[0] = RHS[0];
    for (i = 0;  i < n - 1; i++) {
        k = i + 1;
        l[i] = lower[i] / u[i];
        u[k] = mid[k] - l[i] * upper[i];
        z[k] = RHS[k] - l[i] * z[i];
    }

    ans[n-1] = z[n-1] / u[n-1];
    for (i = 2; i <= n; i++) {
        k = n - i;
        ans[k] = (z[k] - upper[k] * ans[k+1]) / u[k];
    }

    // Free local allocated memory
    free(u);
    free(l);
    free(z);
}

void triCyclicLU(int n, Carray upper, Carray lower, Carray mid, Carray RHS,
     Carray ans)
{
    // Intermediate steps - L . U decomposition
    // Additional line in L denoted by g vector
    // Additional column in U defined bt h vector
    Carray u = carrDef(n);
    Carray l = carrDef(n - 2);
    Carray g = carrDef(n - 1);
    Carray h = carrDef(n - 1);
    Carray z = carrDef(n);
    
    // Loop counters variables.
    unsigned int i, k;
    
    u[0] = mid[0];
    z[0] = RHS[0];
    /****** extras steps ******/
    g[0] = lower[n-1] / mid[0];
    h[0] = upper[n-1];
    /*************************/
    for (i = 0;  i < n - 2; i++) {
        k = i + 1;
        l[i] = lower[i] / u[i];
        u[k] = mid[k] - l[i] * upper[i];
        z[k] = RHS[k] - l[i] * z[i];
        /****** extras steps ******/
        h[k] = (-1) * l[i] * h[i];
        g[k] = (-1) * upper[i] * g[i] / u[k];
        /*************************/
    }
    // little correction due to last terms - L . U = A
    h[n-2] += upper[n-2];
    g[n-2] += lower[n-2] / u[n-2];

    // A little change from simple tridiagonal
    z[n-1] = RHS[n-1] - carrDot2(n-1, g, z);
    u[n-1] = mid[n-1] - carrDot2(n-1, g, h);

    ans[n-1] = z[n-1] / u[n-1];
    // store the last column multiplication in U
    carrScalarMultiply(n - 1, h, ans[n-1], h);
    // Additionaly subtract h compared to tridiagonal
    ans[n-2] = (z[n-2] - h[n-2]) / u[n-2];
    for (i = 3; i <= n; i++) {
        k = n - i;
        ans[k] = (z[k] - h[k] - upper[k] * ans[k+1]) / u[k];
    }

    // Free local allocated memory
    free(u);
    free(l);
    free(z);
    free(h);
    free(g);
}

void triCyclicSM(int n, Carray upper, Carray lower, Carray mid,
                 Carray RHS, Carray ans)
{
    // Cyclic elements
    // bottom == lower[n-1];
    // top    == upper[n-1];

    double complex factor;
    // record data from pointer
    double complex recover1 = mid[0],
                   recoverN = mid[n-1];

    Carray x = carrDef(n);
    Carray w = carrDef(n);
    Carray U = carrDef(n);
    Carray V = carrDef(n);

    carrFill(n, 0, U);
    carrFill(n, 0, V);

    // Choice of gamma factor
    if (cabs(mid[0]) == 0) { factor = upper[0]; mid[0] = -factor; }
    else                   { factor = mid[0];   mid[0] = 0;       }
    
    U[0] = factor;
    V[0] = 1;
    U[n-1] = lower[n-1];
    V[n-1] = upper[n-1] / factor;
    // Adjust last main diagonal element
    mid[n-1] -= upper[n-1] * lower[n-1] / factor;

    #pragma omp parallel sections
    {
        #pragma omp section
        triDiag(n, upper, lower, mid, RHS, x);
        #pragma omp section
        triDiag(n, upper, lower, mid, U, w);
    }

    factor = carrDot2(n, V, x) / (1.0 + carrDot2(n, V, w));
    carrUpdate(n, x, (-1) * factor, w, ans);

    mid[0] = recover1;
    mid[n-1] = recoverN;
    // Free local allocated memory
    free(x);
    free(w);
    free(U);
    free(V);
}


            /*****************************************/
            /* TRIDIAGONAL SYSTEMS WITH REAL ENTRIES */
            /*****************************************/


void rtriDiag(int n, Rarray upper, Rarray lower, Rarray mid, Rarray RHS,
              Rarray ans)
{
    // Intermediate steps - L . U decomposition
    Rarray u = rarrDef(n);
    Rarray l = rarrDef(n - 1);
    Rarray z = rarrDef(n);
    
    // Loop counters variables.
    unsigned int i, k;
    
    if (mid[0] == 0 ) {
        // In this case there is a system reduction
        // where we solve for [x1  x3  x4 ... xn]
        // what is equivalent to adjust the two first 
        // equations and starts counters from 1

        double complex RHS1 = RHS[1] - mid[1]   * RHS[0] / upper[0],
                       RHS2 = RHS[2] - lower[1] * RHS[0] / upper[0];

        // u and z factor initizlization with the changed system
        u[1] = lower[0];
        z[1] = RHS1;
        // One iteration need to be performed outer because
        // the change of values in lower and RHS
        l[1] = 0;
        u[2] = mid[2] - l[1] * upper[1];
        z[2] = RHS2 - l[1] * z[1];
        for (i = 2; i < n - 1; i++) {
            k = i + 1;
            l[i] = lower[i] / u[i];
            u[k] = mid[k] - l[i] * upper[i];
            z[k] = RHS[k] - l[i] * z[i];
        }

        ans[n-1] = z[n-1] / u[n-1];
        for (i = 2; i <= n - 1; i++) {
            k = n - i;
            ans[k] = (z[k] - upper[k] * ans[k+1]) / u[k];
        }

        // Obtained order ans[0..n] = [nan  x1  x3  x4 .. xn]
        // Organize ans[0..n] = [x1  x2  x3  .. xn]
        ans[0] = ans[1];
        ans[1] = RHS[0] / upper[0];

        // Free local alilocated memory
        free(u);
        free(l);
        free(z);

        return;
    }
    
    u[0] = mid[0];
    z[0] = RHS[0];
    for (i = 0;  i < n - 1; i++) {
        k = i + 1;
        l[i] = lower[i] / u[i];
        u[k] = mid[k] - l[i] * upper[i];
        z[k] = RHS[k] - l[i] * z[i];
    }

    ans[n-1] = z[n-1] / u[n-1];
    for (i = 2; i <= n; i++) {
        k = n - i;
        ans[k] = (z[k] - upper[k] * ans[k+1]) / u[k];
    }

    // Free local allocated memory
    free(u);
    free(l);
    free(z);
}
