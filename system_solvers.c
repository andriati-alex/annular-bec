/****** Header file ******/
#include "system_solvers.h"

double cond(unsigned int n, Carray upper, Carray lower, Carray mid)
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
    for (i = 0;  i < n - 1; i++) {
        l = lower[i] / u;
        u = mid[i+1] - l * upper[i];
        condK = 1 + cabs(upper[i] * l / u) * (2 + condK);
        if (maxCond < condK) { maxCond = condK; }
    }

    return maxCond;
}

double errBack(unsigned int n, Carray upper, Carray lower, Carray mid)
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

double complex checkLU(unsigned int n, Carray upper, Carray lower, Carray mid)
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

void triDiag(unsigned int n, Carray upper, Carray lower, Carray mid,
             Carray RHS, Carray ans)
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

        double complex aux = lower[1];

        RHS[1] -= mid[1]   * RHS[0] / upper[0];
        RHS[2] -= lower[1] * RHS[0] / upper[0];

        u[1] = lower[0];
        z[1] = RHS[1];
        lower[1] = 0; // Change value to the new system
        for (i = 1;  i < n - 1; i++) {
            k = i + 1;
            l[i] = lower[i] / u[i];
            u[k] = mid[k] - l[i] * upper[i];
            z[k] = RHS[k] - l[i] * z[i];
        }
        lower[1] = aux; // return original value
        
        // Restore to do not damage outer pointers
        RHS[1] += mid[1]   * RHS[0] / upper[0];
        RHS[2] += lower[1] * RHS[0] / upper[0];

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

void triCyclicLU(unsigned int n, Carray upper, Carray lower, Carray mid,
                 Carray RHS, Carray ans)
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
    carrScalar(n - 1, h, ans[n-1], h);
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

void triCyclicSM(unsigned int n, Carray upper, Carray lower, Carray mid,
                 Carray RHS, Carray ans)
{
    // Cyclic elements
    // bottom == lower[n-1];
    // top    == upper[n-1];

    double complex factor;

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

    triDiag(n, upper, lower, mid, RHS, x);
    triDiag(n, upper, lower, mid, U, w);

    factor = carrDot2(n, V, x) / (1.0 + carrDot2(n, V, w));
    carrScalar(n, w, factor, w);
    carrSub(n, x, w, ans);

    // Free local allocated memory
    free(x);
    free(w);
    free(U);
    free(V);
}

int CCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;
    int l = 0; // Interation counter - return for convergence statistics

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues

    cmatvec(n, n, A, x, r);
    carrSub(n, b, r, r);
    carrCopy(n, r, d);
    while (carrMod(n, r) > eps) {
        cmatvec(n, n, A, d, aux);   // A matrix-vector multiplication
        a = carrMod2(n, r) / carrDot(n, d, aux);
        carrScalar(n, aux, a, aux);
        carrCopy(n, r, prev_r);     // Store (i-1)-th residue to compute beta
        carrSub(n, r, aux, r);      // Update residue
        carrScalar(n, d, a, aux);
        carrAdd(n, x, aux, x);      // Update solution
        beta = carrMod2(n, r) / carrMod2(n, prev_r);
        carrScalar(n, d, beta, aux);
        carrAdd(n, r, aux, d);      // Update direction
        l = l + 1;                  // Update iterarion counter
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_r);

    return l;
}

int preCCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps, 
           Cmatrix M)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;
    int l = 0; // Interation counter - return for convergence statistics

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    cmatvec(n, n, A, x, r);
    carrSub(n, b, r, r);
    cmatvec(n, n, M, r, d);
    carrCopy(n, d, M_r);
    while (carrMod(n, r) > eps) {
        cmatvec(n, n, A, d, aux);   // A matrix-vector multiplication
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        carrScalar(n, aux, a, aux);
        carrCopy(n, r, prev_r);     // Store (i-1)-th residue to compute beta
        carrSub(n, r, aux, r);      // Update residue
        carrScalar(n, d, a, aux);
        carrAdd(n, x, aux, x);      // Update solution
        carrCopy(n, M_r, aux);      // aux get M-1 . r
        cmatvec(n, n, M, r, M_r);   // Store for the next loop
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        carrScalar(n, d, beta, aux);
        carrAdd(n, M_r, aux, d);    // Update direction
        l = l + 1;                  // Update iteration counter
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_r);
    free(M_r);

    return l;
}

int MpreCCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps, 
            Carray upper, Carray lower, Carray mid)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;
    int l = 0; // Interation counter - return for convergence statistics

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    cmatvec(n, n, A, x, r);
    carrSub(n, b, r, r);
    triDiag(n, upper, lower, mid, r, d);
    carrCopy(n, d, M_r);
    while (carrMod(n, r) > eps) {
        cmatvec(n, n, A, d, aux);   // A matrix-vector multiplication
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        carrScalar(n, aux, a, aux);
        carrCopy(n, r, prev_r);     // Store (i-1)-th residue to compute beta
        carrSub(n, r, aux, r);      // Update residue
        carrScalar(n, d, a, aux);
        carrAdd(n, x, aux, x);      // Update solution
        carrCopy(n, M_r, aux);      // aux get M-1 . r
        triDiag(n, upper, lower, mid, r, M_r); // Store for the next loop
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        carrScalar(n, d, beta, aux);
        carrAdd(n, M_r, aux, d);    // Update direction
        l = l + 1;                  // Update iteration counter
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_r);
    free(M_r);

    return l;
}

int CCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;
    int l = 0; // Interation counter - return for convergence statistics

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues

    carrFill(n, 0, r); // CCSvec needs array filled with zeros
    CCSvec(n, A->vec, A->col, A->m, x, r);
    carrSub(n, b, r, r);
    carrCopy(n, r, d);
    while (carrMod(n, r) > eps) {
        carrFill(n, 0, aux); // CCSvec needs array filled with zeros
        CCSvec(n, A->vec, A->col, A->m, d, aux); // matrix-vector mult
        a = carrMod2(n, r) / carrDot(n, d, aux);
        carrScalar(n, aux, a, aux);
        carrCopy(n, r, prev_r);     // Store (i-1)-th residue to compute beta
        carrSub(n, r, aux, r);      // Update residue
        carrScalar(n, d, a, aux);
        carrAdd(n, x, aux, x);      // Update solution
        beta = carrMod2(n, r) / carrMod2(n, prev_r);
        carrScalar(n, d, beta, aux);
        carrAdd(n, r, aux, d);      // Update direction
        l = l + 1;                  // Update iterarion counter
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_r);

    return l;
}

int preCCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps,
              Cmatrix M)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;
    int l = 0; // Interation counter - return for convergence statistics

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    carrFill(n, 0, r); // CCSvec needs array filled with zeros
    CCSvec(n, A->vec, A->col, A->m, x, r);
    carrSub(n, b, r, r);
    cmatvec(n, n, M, r, d);
    carrCopy(n, d, M_r);
    while (carrMod(n, r) > eps) {
        carrFill(n, 0, aux); // CCSvec needs array filled with zeros
        CCSvec(n, A->vec, A->col, A->m, d, aux); // matrix-vector mult
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        carrScalar(n, aux, a, aux);
        carrCopy(n, r, prev_r);     // Store (i-1)-th residue to compute beta
        carrSub(n, r, aux, r);      // Update residue
        carrScalar(n, d, a, aux);
        carrAdd(n, x, aux, x);      // Update solution
        carrCopy(n, M_r, aux);      // aux get M-1 . r
        cmatvec(n, n, M, r, M_r);   // Store for the next loop
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        carrScalar(n, d, beta, aux);
        carrAdd(n, M_r, aux, d);    // Update direction
        l = l + 1;                  // Update iteration counter
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_r);
    free(M_r);

    return l;
}

int MpreCCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps, 
            Carray upper, Carray lower, Carray mid)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;
    int l = 0; // Interation counter - return for convergence statistics

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    carrFill(n, 0, r); // CCSvec needs array filled with zeros
    CCSvec(n, A->vec, A->col, A->m, x, r);
    carrSub(n, b, r, r);
    triDiag(n, upper, lower, mid, r, d);
    carrCopy(n, d, M_r);
    while (carrMod(n, r) > eps) {
        carrFill(n, 0, aux); // CCSvec needs array filled with zeros
        CCSvec(n, A->vec, A->col, A->m, d, aux); // matrix-vector mult
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        carrScalar(n, aux, a, aux);
        carrCopy(n, r, prev_r);     // Store (i-1)-th residue to compute beta
        carrSub(n, r, aux, r);      // Update residue
        carrScalar(n, d, a, aux);
        carrAdd(n, x, aux, x);      // Update solution
        carrCopy(n, M_r, aux);      // aux get M-1 . r
        triDiag(n, upper, lower, mid, r, M_r); // Store for the next loop
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        carrScalar(n, d, beta, aux);
        carrAdd(n, M_r, aux, d);    // Update direction
        l = l + 1;                  // Update iteration counter
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_r);
    free(M_r);

    return l;
}
