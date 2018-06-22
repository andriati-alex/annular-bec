#include "../include/iterative_solver.h"

int CCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;
    int l = 0; // Interation counter - return for convergence statistics

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_x = carrDef(n);
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues

    cmatvec(n, n, A, x, aux);
    carrSub(n, b, aux, r);
    carrCopy(n, r, d);
    while (carrMod(n, r) > eps) {
        // Matrix-Vector multiply
        cmatvec(n, n, A, d, aux);
        a = carrMod2(n, r) / carrDot(n, d, aux);
        // Update residue
        carrCopy(n, r, prev_r);
        carrUpdate(n, prev_r, -a, aux, r);
        // Update solution
        carrCopy(n, x, prev_x);
        carrUpdate(n, prev_x, a, d, x);
        beta = carrMod2(n, r) / carrMod2(n, prev_r);
        // Update direction
        carrScalarMultiply(n, d, beta, aux);
        carrAdd(n, r, aux, d);
        l = l + 1;
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_x);
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
    Carray aux  = carrDef(n);   // to store some algebra
    Carray prev_x = carrDef(n);
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    cmatvec(n, n, A, x, aux);
    carrSub(n, b, aux, r);
    cmatvec(n, n, M, r, d);
    carrCopy(n, d, M_r);
    while (carrMod(n, r) > eps) {
        // Matrix-Vector multiply
        cmatvec(n, n, A, d, aux);
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        // Update residue
        carrCopy(n, r, prev_r);
        carrUpdate(n, prev_r, -a, aux, r);
        // Update solution
        carrCopy(n, x, prev_x);
        carrUpdate(n, prev_x, a, d, x);
        // To compute beta
        carrCopy(n, M_r, aux);
        cmatvec(n, n, M, r, M_r);
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        // Update direction
        carrScalarMultiply(n, d, beta, aux);
        carrAdd(n, M_r, aux, d);
        l = l + 1;
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_x);
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
    Carray prev_x = carrDef(n);
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    cmatvec(n, n, A, x, aux);
    carrSub(n, b, aux, r);
    triDiag(n, upper, lower, mid, r, d);
    carrCopy(n, d, M_r);
    while (carrMod(n, r) > eps) {
        // Matrix-vector multiply
        cmatvec(n, n, A, d, aux);
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        // Update residue
        carrCopy(n, r, prev_r);
        carrUpdate(n, prev_r, -a, aux, r);
        // Update solution
        carrCopy(n, x, prev_x);
        carrUpdate(n, prev_x, a, d, x);
        // to compute beta
        carrCopy(n, M_r, aux);
        triDiag(n, upper, lower, mid, r, M_r);
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        // update direction
        carrScalarMultiply(n, d, beta, aux);
        carrAdd(n, M_r, aux, d);
        l = l + 1;
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_x);
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

    Carray r      = carrDef(n); // residue
    Carray d      = carrDef(n); // Direction
    Carray aux    = carrDef(n); // store intermediate steps
    Carray prev_x = carrDef(n);
    Carray prev_r = carrDef(n);
    
    CCSvec(n, A->vec, A->col, A->m, x, aux);
    carrSub(n, b, aux, r);
    carrCopy(n, r, d);
    
    while (carrMod(n, r) > eps) {
        // Matrix-vector multiply
        CCSvec(n, A->vec, A->col, A->m, d, aux);
        a = carrMod2(n, r) / carrDot(n, d, aux);
        // Update residue
        carrCopy(n, r, prev_r);
        carrUpdate(n, prev_r, -a, aux, r);
        // Update solution
        carrCopy(n, x, prev_x);
        carrUpdate(n, prev_x, a, d, x);
        // Beta and Update direction
        beta = carrMod2(n, r) / carrMod2(n, prev_r);
        carrScalarMultiply(n, d, beta, aux);
        carrAdd(n, r, aux, d);
        l = l + 1;
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_x);
    free(prev_r);

    return l;
}

int preCCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps,
              Cmatrix M)
{
    // Scalars from algorithm
    double complex a;
    double complex beta;

    // Interation counter - return for convergence statistics
    int l = 0;

    Carray r = carrDef(n);      // residue
    Carray d = carrDef(n);      // Direction
    Carray aux = carrDef(n);    // to store some algebra
    Carray prev_x = carrDef(n);
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    CCSvec(n, A->vec, A->col, A->m, x, aux);
    carrSub(n, b, aux, r);
    cmatvec(n, n, M, r, d);
    carrCopy(n, d, M_r);
    while (carrMod(n, r) > eps) {
        // Matrix-vector multiply
        CCSvec(n, A->vec, A->col, A->m, d, aux);
        // scalar
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        // Update Residue
        carrCopy(n, r, prev_r);
        carrUpdate(n, prev_r, -a, aux, r);
        // Update solution
        carrCopy(n, x, prev_x);
        carrUpdate(n, prev_x, a, d, x);
        // Process to compute beta
        carrCopy(n, M_r, aux);
        cmatvec(n, n, M, r, M_r);
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        // Update Direction
        carrScalarMultiply(n, d, beta, aux);
        carrAdd(n, M_r, aux, d);
        l = l + 1;
    }

    // Free function allocated memory
    free(r);
    free(d);
    free(aux);
    free(prev_x);
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
    Carray prev_x = carrDef(n);
    Carray prev_r = carrDef(n); // to compute scalars need 2 residues
    Carray M_r = carrDef(n);    // Pre-conditioner applied to residue

    CCSvec(n, A->vec, A->col, A->m, x, aux);
    carrSub(n, b, aux, r);
    triDiag(n, upper, lower, mid, r, d);
    carrCopy(n, d, M_r);
    
    while (carrMod(n, r) > eps) {
        CCSvec(n, A->vec, A->col, A->m, d, aux); // matrix-vector mult
        a = carrDot(n, r, M_r) / carrDot(n, d, aux);
        // Update residue
        carrCopy(n, r, prev_r);
        carrUpdate(n, prev_r, (-1) * a, aux, r);
        // Update solution
        carrCopy(n, x, prev_x);
        carrUpdate(n, prev_x, a, d, x);
        carrCopy(n, M_r, aux);                  // aux get M-1 . r
        triDiag(n, upper, lower, mid, r, M_r);  // Store for the next loop
        beta = carrDot(n, r, M_r) / carrDot(n, prev_r, aux);
        carrScalarMultiply(n, d, beta, aux);
        carrAdd(n, M_r, aux, d); // Update direction
        l = l + 1;               // Update iteration counter
    }

    // Free function local memory
    free(r);
    free(d);
    free(aux);
    free(prev_x);
    free(prev_r);
    free(M_r);

    return l;
}

int RCG(unsigned int n, RCCSmat A, Rarray b, Rarray x, double eps,
        Rarray upper, Rarray lower, Rarray mid)
{
    // Scalars from algorithm
    double a;
    double beta;
    int l = 0; // Interation counter - return for convergence statistics

    Rarray r = rarrDef(n);      // residue
    Rarray d = rarrDef(n);      // Direction
    Rarray aux = rarrDef(n);    // to store some algebra
    Rarray prev_x = rarrDef(n);
    Rarray prev_r = rarrDef(n); // to compute scalars need 2 residues
    Rarray M_r = rarrDef(n);    // Pre-conditioner applied to residue

    RCCSvec(n, A->vec, A->col, A->m, x, aux);
    rarrSub(n, b, aux, r);
    rtriDiag(n, upper, lower, mid, r, d);
    rarrCopy(n, d, M_r);
    
    while (sqrt(rarrDot(n, r, r)) > eps) {
        RCCSvec(n, A->vec, A->col, A->m, d, aux); // matrix-vector mult
        a = rarrDot(n, r, M_r) / rarrDot(n, d, aux);
        // Update residue
        rarrCopy(n, r, prev_r);
        rarrUpdate(n, prev_r, (-1) * a, aux, r);
        // Update solution
        rarrCopy(n, x, prev_x);
        rarrUpdate(n, prev_x, a, d, x);
        rarrCopy(n, M_r, aux);                  // aux get M-1 . r
        rtriDiag(n, upper, lower, mid, r, M_r); // Store for the next loop
        beta = rarrDot(n, r, M_r) / rarrDot(n, prev_r, aux);
        rarrScalarMultiply(n, d, beta, aux);
        rarrAdd(n, M_r, aux, d); // Update direction
        l = l + 1;               // Update iteration counter
    }

    // Free function local memory
    free(r);
    free(d);
    free(aux);
    free(prev_x);
    free(prev_r);
    free(M_r);

    return l;
}
