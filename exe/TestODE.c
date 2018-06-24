#include "../include/rk4.h"
#include "../include/array_operations.h"

/* ********************************************
 *
 * Apply Runge-Kutta to a linear system of ODEs
 *
 * The script ode.py compare the solution using
 * scipy package an show the accordance of both
 * solutions
 *
 * ********************************************/

int main()
{
    /* DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER */
    omp_set_num_threads(1);
    /******************************************************/

    int n = 3;
    double dt = 0.001;
    int N = (int) 10.0 / dt;
    
    double start, time_used; // show time taken

    Carray upper = carrDef(n - 1);
    Carray lower = carrDef(n - 1);
    Carray mid   = carrDef(n);

    mid[0] = -1.4 * I;
    mid[1] = -0.5 * I;
    mid[2] =  0.2 * I;

    upper[0] = -0.8;
    lower[0] =  0.8;
    upper[1] =  1.5;
    lower[1] = -1.5;

    CCSmat F = triToCCS(n, upper, lower, mid);

    Cmatrix S = cmatDef(N + 1, n);
    // Set initial conditions
    S[0][0] = 0.0;
    S[0][1] = 1.0 / sqrt(2.0);
    S[0][2] = 1.0 / sqrt(2.0);

    RK4(n, N, dt, F, S);

    cmat_txt("test_out/test_ode.dat", N + 1, 1, n, 1, S);

    free(upper);
    free(lower);
    free(mid);
    CCSFree(F);
    cmatFree(N, S);

    printf("\n");
    return 0;
}
