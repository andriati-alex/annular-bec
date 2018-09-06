#include "../include/MCTDHB_integrator"

double complex Energy(MCTDHBsetup mc, Cmatrix Orb, Carray C)
{
    int Morb = mc->Morb;

    int
        i,
        j,
        k,
        l,
        s,
        q;

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
    cmatFree(Morb, Hint);

    return (w + z / 2);
}
