/*

gcc -o test exe/lapack_inversion_test.c -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lm -fopenmp -lmctdhb -L./lib -I./include

*/


#include "../include/MCTDHB_configurations.h"
#include "../include/array_memory.h"


int main(int argc, char * argv[])
{

    omp_set_num_threads(omp_get_max_threads() / 2);
    mkl_set_num_threads(omp_get_max_threads() / 2);


    /* ==================================================================== *
     *                                                                      *
     *                      VARIABLES AND WHAT THEY DO                      *
     *                                                                      *
     * ==================================================================== */

    int
        i,
        j,
        k,
        s,
        Npar, // # of particles
        Morb, // # of orbitals
        ** NCmat,
        ** IF;

    double
        x;

    double complex
        z;

    Carray
        C,
        rho2;

    Cmatrix
        rho;





    /* CASE 1 */

    printf("\n\n==================================================\n\n");
    printf("Case Npar = 2 and Morb = 4.\n\n");

    Npar = 2;
    Morb = 4;

    rho2 = carrDef(Morb * Morb * Morb * Morb);

    rho  = cmatDef(Morb, Morb);

    NCmat = MountNCmat(Npar, Morb);

    IF = MountFocks(Npar, Morb, NCmat);

    C = carrDef(NC(Npar, Morb));

    // Print Configurations and their coefficient index

    printf("\n\nTotal Number of configurations: %d", NC(Npar, Morb));

    printf("\n\nNCmat\n");

    for (i = 0; i < Npar + 1; i++)
    {
        printf("\n\t|");
        for (j = 0; j < Morb + 1; j++) printf(" %5d |", NCmat[i][j]);
    }

    printf("\n\nConfigurations and their indexes\n");
    for (i = 0; i < NCmat[Npar][Morb]; i++)
    {
        printf("\n\t%5d : |", i);
        for (j = 0; j < Morb; j++) printf(" %2d |", IF[i][j]);
        printf(" : %5d", FockToIndex(Npar, Morb, NCmat, IF[i]));
    }

    carrFill(NC(Npar, Morb), 0, C);
    C[3] = 1.0 / sqrt(2);
    C[7] = - 1.0 / sqrt(2) * I;

    OBrho(Npar, Morb, NCmat, IF, C, rho);

    TBrho(Npar, Morb, NCmat, IF, C, rho2);
    
    printf("\n\nOne-Body density matrix\n");
    cmat_print(Morb, Morb, rho);


    printf("\n\nTwo-Body density matrix\n");
    for (i = 0; i < Morb; i ++)
    {
        for (j = 0; j < Morb; j++)
        {
            for (k = 0; k < Morb; k++)
            {
                for (s = 0; s < Morb; s++)
                {
                    z = rho2[i + j*Morb + k*Morb*Morb + s*Morb*Morb*Morb];
                    if (cabs(z) > 1E-10)
                    {
                        printf("\n");
                        printf("rho2[%d, %d, %d, %d] = ", i, j, k, s);
                        cPrint(z);
                    }
                }
            }
        }
    }

    for (i = 0; i < Npar + 1; i++) free(NCmat[i]);
    free(NCmat);
    for (i = 0; i < NC(Npar,Morb); i++) free(IF[i]);
    free(IF);
    cmatFree(Morb, rho);
    free(rho2);
    free(C);
    
    /* CASE 2 */

    printf("\n\n==================================================\n\n");
    printf("Case Npar = 4 and Morb = 3.\n\n");

    Npar = 4;
    Morb = 3;

    rho2 = carrDef(Morb * Morb * Morb * Morb);

    rho  = cmatDef(Morb, Morb);

    NCmat = MountNCmat(Npar, Morb);

    IF = MountFocks(Npar, Morb, NCmat);

    C = carrDef(NC(Npar, Morb));

    // Print Configurations and their coefficient index

    printf("\n\nTotal Number of configurations: %d", NC(Npar, Morb));

    printf("\n\nNCmat\n");

    for (i = 0; i < Npar + 1; i++)
    {
        printf("\n\t|");
        for (j = 0; j < Morb + 1; j++) printf(" %5d |", NCmat[i][j]);
    }

    printf("\n\nConfigurations and their indexes\n");
    for (i = 0; i < NCmat[Npar][Morb]; i++)
    {
        printf("\n\t%5d : |", i);
        for (j = 0; j < Morb; j++) printf(" %2d |", IF[i][j]);
        printf(" : %5d", FockToIndex(Npar, Morb, NCmat, IF[i]));
    }

    carrFill(NC(Npar, Morb), 0, C);
    C[3] = 1.0 / sqrt(2);
    C[7] = - 1.0 / sqrt(2) * I;

    OBrho(Npar, Morb, NCmat, IF, C, rho);

    TBrho(Npar, Morb, NCmat, IF, C, rho2);
    
    printf("\n\nOne-Body density matrix\n");
    cmat_print(Morb, Morb, rho);


    printf("\n\nTwo-Body density matrix\n");
    for (i = 0; i < Morb; i ++)
    {
        for (j = 0; j < Morb; j++)
        {
            for (k = 0; k < Morb; k++)
            {
                for (s = 0; s < Morb; s++)
                {
                    z = rho2[i + j*Morb + k*Morb*Morb + s*Morb*Morb*Morb];
                    if (cabs(z) > 1E-10)
                    {
                        printf("\n");
                        printf("rho2[%d, %d, %d, %d] = ", i, j, k, s);
                        cPrint(z);
                    }
                }
            }
        }
    }
    
    for (i = 0; i < Npar + 1; i++) free(NCmat[i]);
    free(NCmat);
    for (i = 0; i < NC(Npar,Morb); i++) free(IF[i]);
    free(IF);
    cmatFree(Morb, rho);
    free(rho2);
    free(C);

    /* CASE 3 */

    printf("\n\n==================================================\n\n");
    printf("Case Npar = 7 and Morb = 8.\n\n");

    Npar = 7;
    Morb = 8;

    rho2 = carrDef(Morb * Morb * Morb * Morb);

    rho  = cmatDef(Morb, Morb);

    NCmat = MountNCmat(Npar, Morb);

    IF = MountFocks(Npar, Morb, NCmat);

    C = carrDef(NC(Npar, Morb));

    // Print Configurations and their coefficient index

    printf("\n\nTotal Number of configurations: %d", NC(Npar, Morb));

    printf("\n\nNCmat\n");

    for (i = 0; i < Npar + 1; i++)
    {
        printf("\n\t|");
        for (j = 0; j < Morb + 1; j++) printf(" %5d |", NCmat[i][j]);
    }

    printf("\n\nConfigurations and their indexes\n");
    for (i = 0; i < NCmat[Npar][Morb]; i++)
    {
        printf("\n\t%5d : |", i);
        for (j = 0; j < Morb; j++) printf(" %2d |", IF[i][j]);
        j = FockToIndex(Npar, Morb, NCmat, IF[i]);
        printf(" : %5d", j);
        if (j != i)
        {
            printf("\n\n\t\tFATAL ERROR!");
            return -1;
        }
        j = NC(Npar,Morb);
        C[i] = (2 * i - j * I) * cos(PI * ((double) i) / j);
    }

    renormalizeVector(NC(Npar, Morb), C, 1.0);

    OBrho(Npar, Morb, NCmat, IF, C, rho);

    TBrho(Npar, Morb, NCmat, IF, C, rho2);

    carr_txt("test_C.dat", NC(Npar, Morb), C);
    carr_txt("test_rho2.dat", Morb * Morb * Morb * Morb, rho2);
    cmat_txt("rho.dat", Morb, 1, Morb, 1, rho);



    /* ==================================================================== *
     *                                                                      *
     *                      FINISH UP - RELEASE MEMORY                      *
     *                                                                      *
     * ==================================================================== */

    for (i = 0; i < Npar + 1; i++) free(NCmat[i]);
    free(NCmat);
    for (i = 0; i < NC(Npar,Morb); i++) free(IF[i]);
    free(IF);
    cmatFree(Morb, rho);
    free(rho2);
    free(C);

    printf("\n\n");
    return 0;
}
