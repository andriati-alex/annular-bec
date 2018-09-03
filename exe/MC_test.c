#include "../include/MCTDHB_integrator.h"

int main(int argc, char * argv[])
{
    omp_set_num_threads(omp_get_max_threads() / 2);

    int Npar,  // # of particles
        Morb,  // # of orbitals
        ** IF; // Map of index to fock vectors

    int k, l, q, s; // Orbital Counter

    long ** NCmat;  // NCmat[n][m] all combinations of n particles in m boxes
    
    /* ==================================================================== *
     *                                                                      *
     *          SOME VALUES OF ONE- AND TWO-BODY DENSITY MATRICES           *
     *                                                                      *
     * ==================================================================== */

    Npar = 3;
    Morb = 3;
    
    NCmat = MountNCmat(Npar, Morb);

    IF = MountFocks(Npar, Morb, NCmat);

    Carray C = carrDef(NC(Npar, Morb));

    // Setup Values for coeficients. There are 10 altogether
    carrFill(NC(Npar, Morb), 1, C);
    C[2] = - 1 * I; C[5] = 1 + 1 * I; C[6] = 1 - 1 * I; C[8] = 0;

    // One/Two-Body Density Matrices
    Cmatrix rho = cmatDef(Morb, Morb);
    Carray rho2 = carrDef(Morb * Morb * Morb * Morb);



    /*                        CALL ROUTINES                     *
     * -------------------------------------------------------- *
     * Call it twice to see if there is any data race condition */



    OBrho(Npar, Morb, NCmat, IF, C, rho);
    OBrho(Npar, Morb, NCmat, IF, C, rho);

    TBrho(Npar, Morb, NCmat, IF, C, rho2);
    TBrho(Npar, Morb, NCmat, IF, C, rho2);

    printf("\n\n");
    printf("\t=========================================================\n\n");
    printf("\t     Values of Density matrices for Morb = Npar = 3      \n\n");
    printf("\t=========================================================\n\n");

    printf("One-Body: rho \n");
    printf("------------- \n");
    for (k = 0; k < Morb; k++)
    {
        printf("\n\t");
        for (l = 0; l < Morb; l++)
        {
            cPrint(rho[k][l]);
            printf("  ");
        }
    }

    printf("\n\n");
    printf("Two-Body: rho2\n");
    printf("--------------\n");

    for (k = 0; k < Morb; k++)
    {
        for (l = 0; l < Morb; l++)
        {
            for (q = 0; q < Morb; q++)
            {
                for (s = 0; s < Morb; s++)
                {
                    printf("\n\trho2[%d,%d,%d,%d] = ", k, l, q, s);
                    cPrint(rho2[k + l*Morb + q*Morb*Morb + s*Morb*Morb*Morb]);
                }
            }
        }
    }
    
    for (k = 0; k < Npar + 1; k++) free(NCmat[k]); free(NCmat);
    
    for (k = 0; k < NC(Npar, Morb); k++) free(IF[k]); free(IF);

    cmatFree(Morb, rho);

    free(rho2);

    free(C);

    printf("\n\n");
    return 0;
}
