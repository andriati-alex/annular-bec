#include <string.h>
#include "../include/MCTDHB_integrator.h"

void TimePrint(double t)
{
    int
        tt = (int) t,
        days  = 0,
        hours = 0,
        mins  = 0;

    if ( tt / 86400 > 0 )
    { days  = tt / 86400; tt = tt % 86400; }
    if ( tt / 3600  > 0 )
    { hours = tt / 3600;  tt = tt % 3600;  }
    if ( tt / 60 > 0 )
    { mins  = tt / 60;    tt = tt % 60;    }

    printf("%d day(s) %d hour(s) %d minute(s)", days, hours, mins);
}

int main(int argc, char * argv[])
{
    omp_set_num_threads(omp_get_max_threads());


    if (argc != 5)
    {
        printf("\nInvalid Number of command line arguments, ");
        printf("expected: 4.\n\n");
        return -1;
    }

    /* ==================================================================== *
     *                                                                      *
     *                      VARIABLES AND WHAT THEY DO                      *
     *                                                                      *
     * ==================================================================== */

    int // Counters
        k,
        l,
        s,
        what_todo; // real or imaginary time

    int
        N,    // # of time steps to evolve the system
        Mdx,  // # of divisions in space (# of points - 1)
        Npar, // # of particles
        Morb; // # of orbitals

    double
        start, // to measure time
        time_used,
        xi,
        xf,    // Domain of orbitals [xi, xf]
        dt,    // time step (both for real and imaginary)  
        real,  // real part of read data from file
        imag,  // imag part of read data from file
        a2,    // Term multiplying d2 / dx2
        inter, // contact interaction strength
        * V,   // Potential computed in discretized positions
        check; // To check orthogonality condition

    double complex
        a1,        // Term multiplying d / dx
        checkDiag; // To check normalization condition

    char // file configuration name
        fname_in[40],
        fname_out[40];
    
    FILE // pointer to file opened
        * eq_setup_file,
        * out_data;

    Carray
        Ctimetest,
        C,      // Coeficients of superposition of Fock states
        to_int, // auxiliar to compute integration
        E;      // Energy at each time step evolved

    Cmatrix
        Orbtimetest,
        Orb;    // Orb[k][j] give the value of k-th orbital at position j



    /* ==================================================================== *
     *                                                                      *
     *                    READ ALL DATA FROM SETUP FILES                    *
     *                                                                      *
     * ==================================================================== */



    printf("\n\n\n");
    printf("\t=========================================================\n\n");
    printf("\t    Use MC_%s setup files to configure integrator\n\n", argv[3]);
    printf("\t=========================================================\n\n");


    /* Setup Npar, Morb and Mdx *
     * ------------------------ */


    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_config.dat");

    printf("\n\nLooking for %s", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL) // impossible to open file
    { printf("\nERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ... Found !"); }

    k = fscanf(eq_setup_file, "%d %d %d %lf %lf",
               &Npar, &Morb, &Mdx, &xi, &xf);

    fclose(eq_setup_file); // finish reading of file


    /* Setup orbitals *
     * -------------- */


    Orb = cmatDef(Morb, Mdx + 1);
    Orbtimetest = cmatDef(Morb, Mdx + 1);

    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_orb.dat");

    printf("\nLooking for %s ", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("\nERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ... Found !"); }

    for (k = 0; k < Mdx + 1; k++)
    {
        for (s = 0; s < Morb; s++)
        {
            l = fscanf(eq_setup_file, " (%lf+%lfj) ", &real, &imag);
            Orb[s][k] = real + I * imag;
            Orbtimetest[s][k] = real + I * imag;
        }
    }

    fclose(eq_setup_file); // finish the reading of file


    /* Setup Coeficients *
     * ----------------- */


    C = carrDef(NC(Npar, Morb));
    Ctimetest = carrDef(NC(Npar, Morb));
    
    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_coef.dat");

    printf("\nLooking for %s ", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("\nERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ... Found !"); }

    for (k = 0; k < NC(Npar, Morb); k++)
    {
        l = fscanf(eq_setup_file, " (%lf+%lfj)", &real, &imag);
        C[k] = real + I * imag;
        Ctimetest[k] = real + I * imag;
    }

    fclose(eq_setup_file); // finish the reading of file


    /* Setup Equation parameters *
     * ------------------------- */


    a2 = -0.5;
    a1 =  0.0;
    inter = 1.0;
    V = rarrDef(Mdx + 1);
    rarrFill(Mdx + 1, 0, V);


    /* Setup time-step(dt) and number of time-steps to evolve *
     * ------------------------------------------------------ */


    sscanf(argv[1], "%lf", &dt); // First command line argument
    sscanf(argv[2], "%d",  &N);  // Second command line argument



    /* ==================================================================== *
     *                                                                      *
     *                          CALL THE INTEGRATOR                         *
     *                                                                      *
     * ==================================================================== */



    printf("\n\n\n");
    printf("\t=========================================================\n\n");
    printf("\t     Configuration done. Checking orthonormalization     \n\n");
    printf("\t=========================================================\n\n");


    MCTDHBsetup mc = AllocMCTDHBdata(Npar, Morb, Mdx + 1, xi, xf,
                                     a2, inter, V, a1);

    to_int = carrDef(Mdx + 1);


    /* Check if off-diagonal elements are zero *
     * --------------------------------------- */


    check = 0;
    for (k = 0; k < Morb; k++)
    {
        for (l = 0; l < Morb; l++)
        {
            if (l == k) continue;
            for (s = 0; s < Mdx + 1; s++)
            {
                to_int[s] = conj(Orb[k][s]) * Orb[l][s];
            }
            check += cabs(Csimps(Mdx + 1, to_int, mc->dx));
        }
    }

    if (check > 1E-6)
    {
        printf("\n\n\tNot orthogonal orbitals !\n"); return -1;
    }


    /* Check if Diagonal elements sum up to Morb *
     * ----------------------------------------- */


    checkDiag = 0;
    for (k = 0; k < Morb; k++)
    {
        for (s = 0; s < Mdx + 1; s++)
        {
            to_int[s] = conj(Orb[k][s]) * Orb[k][s];
        }
        checkDiag += Csimps(Mdx + 1, to_int, mc->dx);
    }

    if (abs(creal(checkDiag) - Morb) > 1E-6 || cimag(checkDiag) > 1E-6)
    {
        printf("\n\n\tOrbitals are not normalized to 1 !\n"); return -1;
    }


    /* Check if Diagonal elements sum up to Morb *
     * ----------------------------------------- */


    if ( abs(carrMod2(NC(Npar, Morb), C) - 1) > 1E-6 )
    {
        printf("\n\n\tCoeficients norm is not 1 !\n"); return -1;
    }



    /* ==================================================================== *
     *                                                                      *
     *                          CALL THE INTEGRATOR                         *
     *                                                                      *
     * ==================================================================== */



    printf("\n\n\n");
    printf("\t=========================================================\n\n");
    printf("\t      Calling integrator. May take several(hours?)       \n\n");
    printf("\t=========================================================\n\n");
    
    sscanf(argv[4], "%d", &what_todo); // Last command line argument

    E = carrDef(N);

    if (what_todo == 0)
    {   // First estimate time needed based on 1 step
        printf("\n\nDoing real time propagation ...\n");
        
        start = omp_get_wtime();
        MCTDHB_time_evolution(mc, Orbtimetest, Ctimetest, E, dt, 1, 1);
        time_used = (double) (omp_get_wtime() - start);

        printf("\n\nTime to do 1 step: %.1lf seconds\n", time_used);
        printf("\nTotal time estimated: ");
        TimePrint(time_used * N);
        free(Ctimetest);
        cmatFree(Morb, Orbtimetest);

        printf("\n\n");

        // Start Evolution
        MCTDHB_time_evolution(mc, Orb, C, E, dt, N, 1);
    }
    else
    {   // First estimate time needed based on 1 step
        printf("\n\nDoing imaginary time propagation ...\n");
        
        start = omp_get_wtime();
        MCTDHB_itime_evolution(mc, Orbtimetest, Ctimetest, E, dt, 1, 1);
        time_used = (double) (omp_get_wtime() - start);

        printf("\n\nTime to do 1 step: %.1lf seconds\n", time_used);
        printf("\nTotal time estimated: ");
        TimePrint(time_used * N);
        free(Ctimetest);
        cmatFree(Morb, Orbtimetest);

        printf("\n\n");

        // Start Evolution
        MCTDHB_itime_evolution(mc, Orb, C, E, dt, N, 1);
    }



    /* ==================================================================== *
     *                                                                      *
     *                              Record Data                             *
     *                                                                      *
     * ==================================================================== */


    printf("\n\nRecording data ...");

    // Record Orbital Data
    // -------------------

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[3]);
    if (what_todo == 0)
    { strcat(fname_out, "_orb_time.dat");  }
    else
    { strcat(fname_out, "_orb_itime.dat"); }

    cmat_txt(fname_out, Morb, 1, Mdx + 1, 1, Orb);

    // Record Coeficients Data
    // -----------------------

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[3]);
    if (what_todo == 0)
    { strcat(fname_out, "_coef_time.dat");  }
    else
    { strcat(fname_out, "_coef_itime.dat"); }
    carr_txt(fname_out, mc->nc, C);

    // Record Energy Data
    // ------------------

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[3]);
    if (what_todo == 0)
    { strcat(fname_out, "_E_time.dat");  }
    else
    { strcat(fname_out, "_E_itime.dat"); }
    carr_txt(fname_out, N, E);
    
    // Record Parameters Used
    // ----------------------

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[3]);
    if (what_todo == 0)
    { strcat(fname_out, "_config.dat");  }
    else
    { strcat(fname_out, "_iconfig.dat"); }

    out_data = fopen(fname_out, "w");

    if (out_data == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_out); return -1; }

    fprintf(out_data, "%d %d %d %.5lf %.5lf %.10lf %d", 
            Npar, Morb, Mdx, xi, xf, dt, N);

    fclose(out_data);



    /* ==================================================================== *
     *                                                                      *
     *                      FINISH UP - RELEASE MEMORY                      *
     *                                                                      *
     * ==================================================================== */



    EraseMCTDHBdata(mc);
    free(C);
    free(E);
    free(to_int);
    cmatFree(Morb, Orb);

    printf("\n\n");
    return 0;
}
