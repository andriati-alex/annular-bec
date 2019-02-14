#include "data_structure.h"



EqDataPkg PackEqData(int Mpos,double xi,double xf,double a2,double inter,
          doublec a1,char Vname[],double p[])
{

/** Return pointer to a basic data structure with all needed information
  * to solve the Gross-Pitaevskii equation                           **/

    Rarray x = rarrDef(Mpos);

    rarrFillInc(Mpos, xi, (xf - xi) / (Mpos - 1), x);

    EqDataPkg gp = (EqDataPkg) malloc(sizeof(struct _EquationDataPkg));

    if (gp == NULL)
    {
        printf("\n\n\n\tMEMORY ERROR : malloc fail for EqData structure\n\n");
        exit(EXIT_FAILURE);
    }

    gp->Mpos = Mpos;
    gp->xi = xi;
    gp->xf = xf;
    gp->dx = (xf - xi) / (Mpos - 1);
    gp->inter = inter;
    gp->a2 = a2;
    gp->a1 = a1;

    gp->p[0] = p[0];
    gp->p[1] = p[1];
    gp->p[2] = p[2];
    strcpy(gp->Vname,Vname);

    gp->V = rarrDef(Mpos);

    GetPotential(Mpos, Vname, x, gp->V, p[0], p[1], p[2]);

    free(x);

    return gp;
}





void ReleaseEqDataPkg (EqDataPkg gp)
{
    free(gp->V);
    free(gp);
}
