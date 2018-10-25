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





void SetupHint (int Morb, int Mpos, Cmatrix Omat, double dx, double inter,
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
                    Integral = inter * Csimps(Mpos, toInt, dx);
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
