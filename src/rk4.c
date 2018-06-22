#include "../include/rk4.h"

void RK4step(int M, double dt, CCSmat F, Carray y, Carray y_step)
{
    unsigned int i;

    Carray k  = carrDef(M);
    Carray kk = carrDef(M);

    Carray arg = carrDef(M); // arguments to compute k's

    CCSvec(M, F->vec, F->col, F->m, y, k); // k1

    for (i = 0; i < M; i++) {
        arg[i] = 0.5 * dt * k[i] + y[i];  // argument for k2
    }

    CCSvec(M, F->vec, F->col, F->m, arg, kk); // k2

    for (i = 0; i < M; i++) {
        arg[i] = 0.5 * dt * kk[i] + y[i]; // argument for k3
        k[i] += 2 * kk[i];                // Add contribution of  k2
    }

    CCSvec(M, F->vec, F->col, F->m, arg, kk); // k3

    for (i = 0; i < M; i++) {
        arg[i] = dt * kk[i] + y[i]; // argument of k4
        k[i] += 2 * kk[i];          // Add contribution of k3
    }

    CCSvec(M, F->vec, F->col, F->m, arg, kk); // k4

    // Add contribution of k4 and divide by 6 from the method
    for (i = 0; i < M; i++) y_step[i] = y[i] + dt * (k[i] + kk[i]) / 6;

    free(k);
    free(kk);
    free(arg);
}

void RK4(int M, int N, double dt, CCSmat F, Cmatrix S)
{
    unsigned int i, j;
    double r = dt / 6.0;

    Carray k  = carrDef(M);
    Carray kk = carrDef(M);

    Carray arg = carrDef(M); // arguments to compute k's

    for (j = 0; j < N; j++) {

        CCSvec(M, F->vec, F->col, F->m, S[j], k); // k1

        for (i = 0; i < M; i++) {
            arg[i] = 0.5 * dt * k[i] + S[j][i];  // argument for k2
        }

        CCSvec(M, F->vec, F->col, F->m, arg, kk); // k2

        for (i = 0; i < M; i++) {
            arg[i] = 0.5 * dt * kk[i] + S[j][i]; // argument for k3
            k[i] += 2 * kk[i];                   // Add contribution of  k2
        }

        CCSvec(M, F->vec, F->col, F->m, arg, kk); // k3

        for (i = 0; i < M; i++) {
            arg[i] = dt * kk[i] + S[j][i]; // argument of k4
            k[i] += 2 * kk[i];             // Add contribution of k3
        }

        CCSvec(M, F->vec, F->col, F->m, arg, kk); // k4

        // Add contribution of k4 and divide by 6 from the method
        for (i = 0; i < M; i++) S[j+1][i] = S[j][i] + r * (k[i] + kk[i]);
    }

    free(k);
    free(kk);
    free(arg);
}
