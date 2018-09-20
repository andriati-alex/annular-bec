#include "../include/rk4.h"

/* dxdt(size_of_arrays, t, x, extra_argues, derivative_of_x) */

void RK4step(int M, double dt, double t, Carray x, Carray extra, Carray x_step,
     void (*dxdt)(int, double , Carray, Carray, Carray))
{
    int i;

    Carray k = carrDef(M);
    Carray karg = carrDef(M);
    Carray holdk = carrDef(M);

    (*dxdt)(M, t, x, extra, k);

    for (i = 0; i < M; i++)
    {   // initiate with k and prepare argument to compute k2
        holdk[i] = k[i];
        karg[i]  = 0.5 * dt * k[i] + x[i];
    }

    (*dxdt)(M, t + 0.5 * dt, karg, extra, k);

    for (i = 0; i < M; i++)
    {   // Add contribution 2*k2 and prepare arg to compute k3
        holdk[i] += 2 * k[i];
        karg[i]   = 0.5 * dt * k[i] + x[i];
    }

    (*dxdt)(M, t + 0.5 * dt, karg, extra, k);

    for (i = 0; i < M; i++)
    {   // Add contribution 2*k3 and prepare arg to compute k4
        holdk[i] += 2 * k[i];
        karg[i]   = dt * k[i] + x[i];
    }
    
    (*dxdt)(M, t + dt, karg, extra, k);
    
    for (i = 0; i < M; i++)
    {   // Add contribution k4
        holdk[i] += k[i];
    }

    for (i = 0; i < M; i++)
    {   // compute next time step solution
        x_step[i] = x[i] + holdk[i] * dt / 6;
    }

    free(k);
    free(karg);
    free(holdk);
}
