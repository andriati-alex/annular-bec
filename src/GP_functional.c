#include "../include/GP_functional.h"





void applyL0(int n, Carray f, double dx, double a2, double complex a1, 
             Rarray V, double inter, double mu, Carray L0f)
{

/** Apply the nonlinear differential  operator  from  Gross-Pitaevskii
  * stationary equation, that should be zero if function f is solution
  *
  * Output parameter : L0f
  *
**/

    int
        i;


    double complex
        ddx,
        part;


    double
        abs2,
        d2dx;


    ddx  = a1 / (2 * dx);
    d2dx = a2 / (dx * dx);


    #pragma omp parallel for private(i, abs2, part)
    for (i = 1; i < n - 1; i++)
    {

        abs2 = creal(f[i]) * creal(f[i]) + cimag(f[i]) * cimag(f[i]);

        // add first order derivative contribution
        part = (f[i + 1] - f[i - 1]) * ddx;

        // add second order derivative contribution
        part = part + (f[i + 1] - 2 * f[i] + f[i - 1]) * d2dx;

        // add potential contribution
        part = part + f[i] * (V[i] - mu + inter * abs2);

        L0f[i] = part;
    }

    // Finishing using periodic boundary conditions

    abs2 = creal(f[0]) * creal(f[0]) + cimag(f[0]) * cimag(f[0]);
    part =  (f[1] - f[n - 1]) * ddx;
    part = part + (f[1] - 2 * f[0] + f[n - 1]) * d2dx;
    part = part + f[0] * (V[0] - mu + inter * abs2);
    L0f[0] = part;
    
    abs2 = creal(f[n-1]) * creal(f[n-1]) + cimag(f[n-1]) * cimag(f[n-1]);
    part = (f[0] - f[n - 2]) * ddx;
    part = part + (f[0] - 2 * f[n-1] + f[n - 2]) * d2dx;
    part = part + f[n-1] * (V[n-1] - mu + inter * abs2);
    L0f[n-1] = part;

}





double complex Functional(int M, double dx, double a2, double complex a1,
               double inter, Rarray V, Carray f)
{

/** Gross-Pitaesvkii functional that yield both energy and chemical
  * potential depending on what is the coefficient that multiply the
  * nonlinearity of fourth order in integration.
  *
  * For energy : inter = g / 2
  * For chemical potential : inter = g
  *
  * REMIND FOR DIRAC DELTA BARRIER
  *
  * The potential part V must then be integrated by trapezium/rectangle
  * rule. Therefore erase V[i] inside the brackets in for loop defining
  * the integrand. Change the following line that define the energy :
  *
  * E = Csimps(M,Int,dx); for (i = 0; i < M; i++) E = E + dx * V[i] * abs2f[i]
  *
  * This shall work for a numerical implementation like 1 / dx
  *
**/

    int
        i;

    double
        norm;

    double complex
        E;

    Carray
        DF  = carrDef(M),
        Int = carrDef(M);

    Rarray
        abs2F  = rarrDef(M),
        abs2DF = rarrDef(M);



    dxFD(M,f,dx,DF);

    carrAbs2(M,DF,abs2DF);

    carrAbs2(M,f,abs2F);

    for (i = 0; i < M; i++)
    {
        // Use E to sum potential and derivatives contributions
        E = (V[i] + inter * abs2F[i]) * abs2F[i];
        E = E - a2 * abs2DF[i] + a1 * conj(f[i]) * DF[i];
        Int[i] = E;
    }


    E = Csimps(M,Int,dx);
    norm = Rsimps(M,abs2F,dx);

    // release memory
    free(DF); free(abs2F); free(abs2DF); free(Int);

    return E / norm;

}





double complex GPkinect(int M, double a2, double complex a1, double dx,
               Carray psi)
{

    int k;

    double complex
        r;

    Carray
        ddx   = carrDef(M),
        toInt = carrDef(M);

    dxFD(M,psi,dx,ddx);

    for (k = 0; k < M; k++)
    {
        toInt[k] =  - a2 * conj(ddx[k]) * ddx[k];
    }

    r = Csimps(M,toInt,dx);

    free(ddx); free(toInt);

    return r;

}





double complex GPtrap(int M, Rarray V, double dx, Carray psi)
{

    int k;

    double complex
        r;

    Carray
        toInt = carrDef(M);

    for (k = 0; k < M; k++)
    {
        toInt[k] = V[k] * conj(psi[k]) * psi[k];
    }

    r = Csimps(M,toInt,dx);

    free(toInt);

    return r;

}





double GPinter(int M, double g, double dx, Carray psi)
{

    int k;

    double
        r;

    Rarray
        toInt = rarrDef(M);

    carrAbs2(M,psi,toInt);

    for (k = 0; k < M; k++)  toInt[k] = g * toInt[k] * toInt[k];

    r = Rsimps(M,toInt,dx) / 2;

    free(toInt);

    return r;

}





double complex GPvirial(int M, double a2, double complex a1, double g,
               Rarray V, double dx, Carray psi)
{
    double complex
        kin,
        pot,
        ntr;

    kin = GPkinect(M,a2,a1,dx,psi);
    pot = GPtrap(M,V,dx,psi);
    ntr = GPinter(M,g,dx,psi);

    return ( 2 * pot - 2 * kin - ntr );
}
