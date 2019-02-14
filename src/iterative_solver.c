#include "iterative_solver.h"





int CCG(int n, CCSmat A, Carray b, Carray x, double eps, int maxiter,
        Carray upper, Carray lower, Carray mid)
{

    double complex
        a,
        beta;

    int
        l; // Interation counter - return for convergence statistics

    Carray
        r,
        d,
        aux,
        prev_x,
        prev_r,
        M_r;



    l = 0;

    r = carrDef(n);      // residue
    d = carrDef(n);      // Direction
    aux = carrDef(n);    // to store some algebra
    prev_x = carrDef(n);
    prev_r = carrDef(n); // to compute scalars need 2 residues
    M_r = carrDef(n);    // Pre-conditioner applied to residue

    CCSvec(n, A->vec, A->col, A->m, x, aux);
    carrSub(n, b, aux, r);
    triDiag(n, upper, lower, mid, r, d);
    carrCopy(n, d, M_r);
    
    while (carrMod(n, r) > eps)
    {
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

        if (l == maxiter)
        {
            printf("\n\n\tWARNING : exit before achieve desired residual ");
            printf("value in Conjugate Gradient method due to max number ");
            printf("of iterations given =  %d\n\n", maxiter);
            break;
        }
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





int RCG(int n, RCCSmat A, Rarray b, Rarray x, double eps, int maxiter,
        Rarray upper, Rarray lower, Rarray mid)
{
    
    double
        a,
        beta;

    int l = 0; // Interation counter - return for convergence statistics

    Rarray
        r,
        d,
        aux,
        prev_x,
        prev_r,
        M_r;

    r = rarrDef(n);      // residue
    d = rarrDef(n);      // Direction
    aux = rarrDef(n);    // to store some algebra
    prev_x = rarrDef(n);
    prev_r = rarrDef(n); // to compute scalars need 2 residues
    M_r = rarrDef(n);    // Pre-conditioner applied to residue

    RCCSvec(n, A->vec, A->col, A->m, x, aux);
    rarrSub(n, b, aux, r);
    realtri(n, upper, lower, mid, r, d);
    rarrCopy(n, d, M_r);
    
    while (sqrt(rarrDot(n, r, r)) > eps)
    {
        RCCSvec(n, A->vec, A->col, A->m, d, aux); // matrix-vector mult
        a = rarrDot(n, r, M_r) / rarrDot(n, d, aux);
        // Update residue
        rarrCopy(n, r, prev_r);
        rarrUpdate(n, prev_r, (-1) * a, aux, r);
        // Update solution
        rarrCopy(n, x, prev_x);
        rarrUpdate(n, prev_x, a, d, x);
        rarrCopy(n, M_r, aux);                  // aux get M-1 . r
        realtri(n, upper, lower, mid, r, M_r); // Store for the next loop
        beta = rarrDot(n, r, M_r) / rarrDot(n, prev_r, aux);
        rarrScalarMultiply(n, d, beta, aux);
        rarrAdd(n, M_r, aux, d); // Update direction

        l = l + 1;               // Update iteration counter

        if (l == maxiter)
        {
            printf("\n\n\tWARNING : exit before achieve desired residual ");
            printf("value in Conjugate Gradient method due to max number ");
            printf("of iterations given =  %d\n\n", maxiter);
            break;
        }
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
