#include "observables.h"





doublec Chem(int M, double dx, double a2, doublec a1, double inter,
        Rarray V, Carray f)
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





doublec KinectE(int M, double a2, doublec a1, double dx, Carray psi)
{

    int k;

    double complex
        r;

    Carray
        ddx   = carrDef(M),
        toInt = carrDef(M);

    dxFD(M,psi,dx,ddx);

    for (k = 0; k < M; k++) toInt[k] = - a2 * conj(ddx[k]) * ddx[k];

    r = Csimps(M,toInt,dx);

    free(ddx); free(toInt);

    return r;

}





doublec TrapE(int M, Rarray V, double dx, Carray psi)
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





double InterE(int M, double g, double dx, Carray psi)
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





doublec Energy(int M, double dx, double a2, doublec a1, double inter,
        Rarray V, Carray f)
{
    return Chem(M, dx, a2, a1, inter / 2, V, f);
}





doublec Virial(int M, double a2, doublec a1, double g, Rarray V,
        double dx, Carray psi)
{
    double complex
        kin,
        pot,
        ntr;

    kin = KinectE(M,a2,a1,dx,psi);
    pot = TrapE(M,V,dx,psi);
    ntr = InterE(M,g,dx,psi);

    return ( 2 * pot - 2 * kin - ntr );
}





double MeanQuadraticR(int n, Carray f, double dx)
{

/** Compute Mean Square value of normalized complex function/distribution **/

    int
        i;

    double
        r;

    Rarray
        ToInt;

    ToInt = rarrDef(n);

    r = - dx * (n - 1) * 0.5; // First discretized point of domain

    for (i = 0; i < n; i ++)
    {
        ToInt[i] = (creal(f[i])*creal(f[i]) + cimag(f[i])*cimag(f[i])) * r*r;
        r = r + dx;
    }

    r = sqrt(Rsimps(n, ToInt, dx));

    free(ToInt);

    return r;

}
