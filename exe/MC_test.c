#include <string.h>
#include "../include/MCTDHB_integrator.h"

int main(int argc, char * argv[])
{
    omp_set_num_threads(omp_get_max_threads());

    int N,     // # of time steps to go foward
        Mdx,   // # of divisions in space
        Npar,  // # of particles
        Morb,  // # of orbitals
        ** IF; // Map of index to fock vectors

    int k, l, q, s; // Orbital Counter

    long ** NCmat;  // NCmat[n][m] all combinations of n particles in m boxes

    double dt, real, imag;



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
    
    
    
    /* ==================================================================== *
     *                                                                      *
     *                    READ ALL DATA FROM SETUP FILES                    *
     *                                                                      *
     * ==================================================================== */



    printf("\n\n\n");
    printf("\t=========================================================\n\n");
    printf("\t     Use MC_%s setup files to call integrator\n\n", argv[3]);
    printf("\t=========================================================\n\n");


    /* Setup Npar, Morb and Mdx *
     * ------------------------ */


    char fname_in[40];    // file configuration names
    FILE * eq_setup_file; // pointer to file

    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_config.dat");

    printf("\n\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL) // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    k = fscanf(eq_setup_file, "%d %d %d", &Npar, &Morb, &Mdx);

    fclose(eq_setup_file); // finish reading of file


    /* Setup time-step(dt) and number of time-steps to evolve *
     * ------------------------------------------------------ */


    sscanf(argv[1], "%lf", &dt); // First command line argument
    sscanf(argv[2], "%d",  &N);  // Second command line argument


    /* Setup orbitals *
     * -------------- */


    Cmatrix Orb = cmatDef(Morb, Mdx + 1);

    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_orb.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    for (k = 0; k < Mdx + 1; k++)
    {
        for (s = 0; s < Morb; s++)
        {
            l = fscanf(eq_setup_file, " (%lf+%lfj) ", &real, &imag);
            Orb[s][k] = real + I * imag;
        }
    }

    fclose(eq_setup_file); // finish the reading of file
    
    /* Setup Coeficients *
     * ----------------- */

    C = carrDef(NC(Npar, Morb));
    
    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_coef.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    for (k = 0; k < NC(Npar, Morb); k++)
    {
        l = fscanf(eq_setup_file, " (%lf+%lfj)", &real, &imag);
        C[k] = real + I * imag;
    }

    fclose(eq_setup_file); // finish the reading of file

    
    
    /* ==================================================================== *
     *                                                                      *
     *                          CALL THE INTEGRATOR                         *
     *                                                                      *
     * ==================================================================== */

    

    double a2    = -0.5,
           inter = 1.0,
           * V   = rarrDef(Mdx + 1);

    rarrFill(Mdx + 1, 0, V);
    
    double complex a1 = 0.0;

    MCTDHBsetup mc = AllocMCTDHBdata(Npar, Morb, Mdx + 1, -PI, PI,
                                     a2, inter, V, a1);

    MCTDHB_time_evolution(mc, Orb, C, dt, N, 1);
    
    /* ==================================================================== *
     *                                                                      *
     *                              Record Data                             *
     *                                                                      *
     * ==================================================================== */
    
    char fname_out[30];

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[3]);
    strcat(fname_out, "_itime.dat");

    printf("\nRecording data ...\n");
    cmat_txt(fname_out, Morb, 1, Mdx + 1, 1, Orb);

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[3]);
    strcat(fname_out, "_config.dat");

    FILE * out_data = fopen(fname_out, "w");

    if (out_data == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_out); return -1; }

    fprintf(out_data, "%d %d %d %.10lf %d", Npar, Morb, Mdx, dt, N);

    fclose(out_data);

    EraseMCTDHBdata(mc);
    free(C);
    cmatFree(Morb, Orb);

    printf("\n\n");
    return 0;
}
