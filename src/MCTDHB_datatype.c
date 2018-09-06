#include "../include/MCTDHB_datatype.h"

MCTDHBsetup AllocMCTDHBdata (
    int Npar, 
    int Morb,
    int Mpos,
    double xi,
    double xf,
    double a2,
    double inter,
    double * V,
    double complex a1 )
{   // Configure and return the pointer to MCTDHB structure

    MCTDHBsetup MC = (MCTDHBsetup) malloc(sizeof(struct _MCTDHBsetup));
    MC->Npar = Npar;
    MC->Morb = Morb;
    MC->Mpos = Mpos;
    MC->xi = xi;
    MC->xf = xf;
    MC->dx = (xf - xi) / (Mpos - 1);
    MC->inter = inter;
    MC->a2 = a2;
    MC->a1 = a1;
    MC->V = V;
    MC->nc = NC(Npar, Morb);
    MC->NCmat = MountNCmat(Npar, Morb);
    MC->IF = MountFocks(Npar, Morb, MC->NCmat);
    return MC;
}

void EraseMCTDHBdata (MCTDHBsetup MC)
{   // release all fields in the structure
    long i;
    for (i = 0; i < MC->nc; i++) free(MC->IF[i]);
    free(MC->IF);
    for (i = 0; i <= MC->Npar; i++) free(MC->NCmat[i]);
    free(MC->NCmat);
    free(MC->V);
    free(MC);
}
