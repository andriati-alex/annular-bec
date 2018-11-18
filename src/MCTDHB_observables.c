#include "../include/MCTDHB_observables.h"





void SetupHo (int Morb, int Mpos, Cmatrix Omat, double dx, double a2,
     double complex a1, Rarray V, Cmatrix Ho )
{   // Setup matrix elements of onebody hamiltonian part

    int i,
        j,
        k;

    Carray ddxi  = carrDef(Mpos);
    Carray ddxj  = carrDef(Mpos);
    Carray toInt = carrDef(Mpos);

    for (i = 0; i < Morb; i++)
    {
        dxCyclic(Mpos, Omat[i], dx, ddxi);
        for (j = 0; j < Morb; j++)
        {
            dxCyclic(Mpos, Omat[j], dx, ddxj);
            for (k = 0; k < Mpos; k++)
            {
                toInt[k] = - a2 * conj(ddxi[k]) * ddxj[k]    \
                           + a1 * conj(Omat[i][k]) * ddxj[k] \
                           + V[k] * conj(Omat[i][k]) * Omat[j][k];
            }
            Ho[i][j] = Csimps(Mpos, toInt, dx);
        }
    }

    free(ddxi); free(ddxj); free(toInt);
}





void SetupHint (int Morb, int Mpos, Cmatrix Omat, double dx, double g,
     Carray Hint)
{   // Matrix elements of twobody hamiltonian part

    int i,
        k,
        s,
        q,
        l,
        M = Morb,
        M2 = Morb * Morb,
        M3 = Morb * Morb * Morb;

    double complex Integral;

    Carray toInt = carrDef(Mpos);

    for (k = 0; k < Morb; k++)
    {
        for (s = k; s < Morb; s++)
        {
            for (q = 0; q < Morb; q++)
            {
                for (l = q; l < Morb; l++)
                {
                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(Omat[k][i] * Omat[s][i]) * \
                                   Omat[q][i] * Omat[l][i];
                    }
                    Integral = g * Csimps(Mpos, toInt, dx);
                    Hint[k + s * M + q * M2 + l * M3] = Integral;
                    Hint[k + s * M + l * M2 + q * M3] = Integral;
                    Hint[s + k * M + q * M2 + l * M3] = Integral;
                    Hint[s + k * M + l * M2 + q * M3] = Integral;
                }   // Take advantage of the symmetry k <--> s
            }
        }
    }

    free(toInt);
}





double complex Energy (MCTDHBsetup mc, Cmatrix Orb, Carray C)
{   // return the energy from multiconfigurational description  with
    // configurations coeficients C respect to occupation in the set
    // of orbitals Orb[k].

    int
        i,
        j,
        k,
        l,
        s,
        q,
        Morb;

    Morb = mc->Morb;

    double complex
        z,
        w;

    Cmatrix
        rho = cmatDef(Morb, Morb),
        Ho  = cmatDef(Morb, Morb);

    Carray
        rho2 = carrDef(Morb * Morb * Morb * Morb),
        Hint = carrDef(Morb * Morb * Morb * Morb);

    OBrho(mc->Npar, Morb, mc->NCmat, mc->IF, C, rho);
    TBrho(mc->Npar, Morb, mc->NCmat, mc->IF, C, rho2);
    
    SetupHo(Morb, mc->Mpos, Orb, mc->dx, mc->a2, mc->a1, mc->V, Ho);
    SetupHint(Morb, mc->Mpos, Orb, mc->dx, mc->inter, Hint);

    z = 0;
    w = 0;
    for (k = 0; k < Morb; k++)
    {
        for (l = 0; l < Morb; l++)
        {
            w += rho[k][l] * Ho[k][l];
            for (s = 0; s < Morb; s++)
            {
                for (q = 0; q < Morb; q++)
                {
                    i = k + l * Morb + s * Morb*Morb + q * Morb*Morb*Morb;
                    j = k + l * Morb + q * Morb*Morb + s * Morb*Morb*Morb;
                    z += rho2[i] * Hint[j];
                }
            }
        }
    }

    free(rho2);
    free(Hint);
    cmatFree(Morb, rho);
    cmatFree(Morb, Ho);

    return (w + z / 2);
}





double complex KinectE (int Morb, int Mpos, Cmatrix Omat, double dx, double a2,
     double complex a1, Cmatrix rho )
{

    int i,
        j,
        k;

    double complex
        r;

    Carray
        ddxi  = carrDef(Mpos),
        ddxj  = carrDef(Mpos),
        toInt = carrDef(Mpos);



    carrFill(Mpos, 0, toInt);
    for (i = 0; i < Morb; i++)
    {
        dxCyclic(Mpos, Omat[i], dx, ddxi);
        for (j = 0; j < Morb; j++)
        {
            r = rho[i][j];
            dxCyclic(Mpos, Omat[j], dx, ddxj);
            for (k = 0; k < Mpos; k++)
            {
                toInt[k] = toInt[k] - a2 * r * conj(ddxi[k]) * ddxj[k];
                toInt[k] = toInt[k] + a1 * r * conj(Omat[i][k]) * ddxj[k];
            }
        }
    }

    r = Csimps(Mpos, toInt, dx);

    free(ddxi); free(ddxj); free(toInt);

    return r;
}





double complex PotentialE (int Morb, int Mpos, Cmatrix Omat, double dx,
     Rarray V, Cmatrix rho )
{

    int i,
        j,
        k;

    double complex
        r;

    Carray
        toInt = carrDef(Mpos);



    carrFill(Mpos, 0, toInt);
    for (i = 0; i < Morb; i++)
    {
        for (j = 0; j < Morb; j++)
        {
            r = rho[i][j];
            for (k = 0; k < Mpos; k++)
                toInt[k] = toInt[k] + r * V[k] * conj(Omat[i][k]) * Omat[j][k];
        }
    }

    r = Csimps(Mpos, toInt, dx);

    free(toInt);

    return r;
}





double complex InteractingE(int Morb, int Mpos, Cmatrix Omat, double dx, double g,
     Carray rho)
{

    int i,
        k,
        s,
        q,
        l,
        M = Morb,
        M2 = Morb * Morb,
        M3 = Morb * Morb * Morb;

    double complex
        r;

    Carray
        toInt = carrDef(Mpos);



    carrFill(Mpos, 0, toInt);
    for (k = 0; k < Morb; k++)
    {

        for (s = 0; s < Morb; s++)
        {

            for (q = 0; q < Morb; q++)
            {

                for (l = 0; l < Morb; l++)
                {

                    r = rho[k + s * M + q * M2 + l * M3];
                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] += r * conj(Omat[k][i] * Omat[s][i]) * \
                        Omat[l][i] * Omat[q][i];
                    }

                }   // Take advantage of the symmetry k/s and q/l
            }
        }
    }

    r = g * Csimps(Mpos, toInt, dx) / 2;

    free(toInt);

    return r;
}





double complex VirialResidue(MCTDHBsetup mc, Cmatrix Orb, Carray C)
{

    int
        Npar = mc->Npar,
        Morb = mc->Morb,
        Mpos = mc->Mpos;

    double
        dx = mc->dx,
        g  = mc->inter,
        a2 = mc->a2,
        *V = mc->V;

    double complex
        kinect,
        potential,
        interacting,
        a1 = mc->a1;
    
    Cmatrix
        rho = cmatDef(Morb, Morb);

    Carray
        rho2 = carrDef(Morb * Morb * Morb * Morb);



    OBrho(Npar, Morb, mc->NCmat, mc->IF, C, rho);
    TBrho(Npar, Morb, mc->NCmat, mc->IF, C, rho2);



    kinect = KinectE(Morb, Mpos, Orb, dx, a2, a1, rho);
    potential = PotentialE(Morb, Mpos, Orb, dx, V, rho);
    interacting = InteractingE(Morb, Mpos, Orb, dx, g, rho2);

    free(rho2);
    cmatFree(Morb, rho);

    return (2 * potential - 2 * kinect - interacting);

}
