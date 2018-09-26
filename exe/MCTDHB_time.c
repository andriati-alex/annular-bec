#include <string.h>
#include "../include/MCTDHB_integrator.h"





/* =========================================================================
 *
 *
 * OBTAIN REAL OR IMAGINARY TIME EVOLUTION SOLUTION OF
 * MCTDHB EQUATIONS WITH DIRAC DELTA BARRIER
 *
 *
 *
 *
 *
 *
 * REQUIRED FILES
 * -------------------------------------------------------------------------
 *
 * setup/MC_fileId_orb.dat
 *
 *      Text file with a matrix where the k-th column represent the k-th
 *      orbital. Thus the rows represent the values of these orbitals in
 *      discretized positions
 *
 * setup/MC_fileId_eq.dat
 *
 *      A text file within the values of equation coefficients in the
 *      following order, separeted by spaces:
 *
 *      (1) second order derivative
 *      (2) imag part of first order derivative(have no real part)
 *      (3) interaction strength (g)
 *      (4) Delta Barrier Strength
 *      (5) Boundary conditions - Boolean
 *
 * setup/MC_fileId_domain.dat
 *
 *      A text file with position domain information over  which  the
 *      fileId_init was generated. The numbers are in following order
 *
 *      (1) # of particles
 *      (2) # of orbitals
 *      (3) x_i
 *      (4) x_f
 *      (5) M the number of slices of size (xf - xi) / M
 *
 *
 *
 *
 *
 * COMMAND LINE ARGUMENTS
 * -------------------------------------------------------------------------
 * 
 * real/imag dt N fileId output_name method(optional)
 *
 *      real/imag   -> real/Real/imag/Imag (the time)
 *      dt          -> the time step
 *      N           -> number of time steps to propagate
 *      fileId      -> Name to look up for files
 *      output_name -> name of file to write results
 *      method      -> Method to solve (see below)
 *
 * PS: method, if given, must have as :
 *
 *      * First algarism define how to do linear part
 *          1 Crank-Nicolson with Sherman-Morrison
 *          2 Crank-Nicolson with LU-decomposition
 *
 *      * Second algarism define how to do nonlinear part
 *          1 4th order Runge-Kutta
 *          2 Mix between 4th order Runge-Kutta and lanczos
 *
 *
 *
 *
 *
 * CALL
 * -------------------------------------------------------------------------
 *
 * ./MCTDHB_time real/imag dt N fileId output_name method(optional)
 *
 *
 *
 *
 *
 * OUTPUT FILES
 * -------------------------------------------------------------------------
 *
 *
 * ========================================================================= */





void TimePrint(double t)
{   // format and print time in days / hours / minutes
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
    omp_set_num_threads(omp_get_max_threads() / 2);


    if (argc < 6 || argc > 7)
    {
        printf("\nInvalid Number of command line arguments, ");
        printf("expected: 5 or 6.\n\n");
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
        s;

    int
        N,    // # of time steps to evolve the system
        Mdx,  // # of divisions in space (# of points - 1)
        Npar, // # of particles
        Morb, // # of orbitals
        cyclic,
        method;

    double
        start,  // to measure time
        time_used,
        dx,
        xi,
        xf,     // Domain of orbitals [xi, xf]
        dt,     // time step (both for real and imaginary)  
        real,   // real part of read data from file
        imag,   // imag part of read data from file
        a2,     // Term multiplying d2 / dx2
        inter,  // contact interaction strength
        lambda, // height of dirac delta barrier
        check,
        * V;   // Potential computed in discretized positions

    double complex
        a1,        // Term multiplying d / dx
        checkDiag; // To check normalization condition

    char // file configuration name
        timeinfo,
        fname_in[120],
        fname_out[120];
    
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
    printf("\t    Use MC_%s setup files to configure integrator\n\n", argv[4]);
    printf("\t=========================================================\n\n");



    /* Setup Npar, Morb and Mdx *
     * ------------------------ */



    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[4]);
    strcat(fname_in, "_config.dat");

    printf("\n\nLooking for %s", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL) // impossible to open file
    { printf("\nERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ... Found !"); }

    k = fscanf(eq_setup_file, "%d %d %d %lf %lf",
               &Npar, &Morb, &Mdx, &xi, &xf);

    dx = (xf - xi) / Mdx;

    printf("\n\t# of Particles: %d", Npar);
    printf("\n\t# of Orbitals: %d", Morb);
    printf("\n\t# of possible configurations: %d\n", NC(Npar, Morb));

    fclose(eq_setup_file); // finish reading of file



    /* Setup orbitals *
     * -------------- */



    Orb = cmatDef(Morb, Mdx + 1);
    Orbtimetest = cmatDef(Morb, Mdx + 1);

    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[4]);
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
    strcat(fname_in, argv[4]);
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



    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[4]);
    strcat(fname_in, "_eq.dat");
    
    printf("\nLooking for %s", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ... Found !"); }    

    l = fscanf(eq_setup_file, "%lf %lf %lf %lf %d",
                   &a2, &imag, &inter, &lambda, &cyclic);

    fclose(eq_setup_file);

    a1 = 0 + imag * I;
    
    V = rarrDef(Mdx + 1);
    rarrFill(Mdx + 1, 0, V);
    V[Mdx/2] = lambda / dx;  // Dirac Delta barrier case lambda != 0

    inter =  - sqrt(2.0) * 4.0 / (Npar - 1);


    /* Setup time-step(dt) and number of time-steps to evolve *
     * ------------------------------------------------------ */



    timeinfo = argv[1][0]; // character i for imaginary or r for real time

    if (timeinfo != 'r' && timeinfo != 'R')
    {
        if (timeinfo != 'i' && timeinfo != 'I')
        {
            printf("\n\n\tInvalid first argument.\n");
            return -1;
        }
    }

    sscanf(argv[2], "%lf", &dt); // First command line argument
    sscanf(argv[3], "%d",  &N);  // Second command line argument
    
    
    
    /* Integrator method *
     * ----------------- */
    if (argc == 7) { sscanf(argv[6], "%d", &method); }
    else           { method = 11;                    }



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

    if (check > 1E-9)
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

    if (abs(creal(checkDiag) - Morb) > 1E-9 || cimag(checkDiag) > 1E-9)
    {
        printf("\n\n\tOrbitals are not normalized to 1 !\n"); return -1;
    }


    /* Check normalization of coeficients *
     * ---------------------------------- */


    if ( abs(carrMod2(NC(Npar, Morb), C) - 1) > 1E-9 )
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

    E = carrDef(N); // Store energy

    if (timeinfo == 'r' || timeinfo == 'R')
    {   // First estimate time needed based on 1 step

        switch (method)
        {
            case 11:
                printf("\n\nDoing real time propagation (SM/RK4) ...\n");
        
                start = omp_get_wtime();
                MCTDHB_REAL_RK4I(mc, Orbtimetest, Ctimetest, E, dt, 1, 1);
                time_used = (double) (omp_get_wtime() - start);

                printf("\n\nTime to do 1 step: %.1lf seconds\n", time_used);
                printf("\nTotal time estimated: ");
                TimePrint(time_used * N);
                free(Ctimetest);
                cmatFree(Morb, Orbtimetest);

                printf("\n\n");

                // Start Evolution
                MCTDHB_REAL_RK4I(mc, Orb, C, E, dt, N, 1);
                break;
            case 12:
                printf("\n\nDoing real time propagation");
                printf(" (SM/RK4/lanczos) ...\n");
        
                start = omp_get_wtime();
                MCTDHB_REAL_LanczosRK4I(mc, Orbtimetest, Ctimetest, E, dt, 1, 1);
                time_used = (double) (omp_get_wtime() - start);

                printf("\n\nTime to do 1 step: %.1lf seconds\n", time_used);
                printf("\nTotal time estimated: ");
                TimePrint(time_used * N);
                free(Ctimetest);
                cmatFree(Morb, Orbtimetest);

                printf("\n\n");

                // Start Evolution
                MCTDHB_REAL_LanczosRK4I(mc, Orb, C, E, dt, N, 1);
                break;
            case 21:
                printf("\n\nDoing real time propagation (LU/RK4) ...\n");
        
                start = omp_get_wtime();
                MCTDHB_REAL_RK4I(mc, Orbtimetest, Ctimetest, E, dt, 1, 1);
                time_used = (double) (omp_get_wtime() - start);

                printf("\n\nTime to do 1 step: %.1lf seconds\n", time_used);
                printf("\nTotal time estimated: ");
                TimePrint(time_used * N);
                free(Ctimetest);
                cmatFree(Morb, Orbtimetest);

                printf("\n\n");

                // Start Evolution
                MCTDHB_REAL_RK4I(mc, Orb, C, E, dt, N, 1);
                break;
            case 22:
                printf("\n\nDoing real time propagation");
                printf(" (LU/RK4/lanczos) ...\n");
        
                start = omp_get_wtime();
                MCTDHB_REAL_LanczosRK4I(mc, Orbtimetest, Ctimetest, E, dt, 1, 1);
                time_used = (double) (omp_get_wtime() - start);

                printf("\n\nTime to do 1 step: %.1lf seconds\n", time_used);
                printf("\nTotal time estimated: ");
                TimePrint(time_used * N);
                free(Ctimetest);
                cmatFree(Morb, Orbtimetest);

                printf("\n\n");

                // Start Evolution
                MCTDHB_REAL_LanczosRK4I(mc, Orb, C, E, dt, N, 1);
                break;
        }
    }
    else
    {
        printf("\n\nDoing imaginary time propagation ...\n");

        start = omp_get_wtime();
        MCTDHB_IMAG_RK4I(mc, Orbtimetest, Ctimetest, E, dt, 1, 1);
        time_used = (double) (omp_get_wtime() - start);

        printf("\n\nTime to do 1 step: %.1lf seconds\n", time_used);
        printf("\nTotal time estimated: ");
        TimePrint(time_used * N);
        free(Ctimetest);
        cmatFree(Morb, Orbtimetest);

        printf("\n\n");

        // Start Evolution
        MCTDHB_IMAG_RK4I(mc, Orb, C, E, dt, N, 1);
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
    strcat(fname_out, argv[5]);
    if (timeinfo == 'r' || timeinfo == 'R')
    { strcat(fname_out, "_orb_realtime.dat");  }
    else
    { strcat(fname_out, "_orb_imagtime.dat"); }

    cmat_txt(fname_out, Morb, 1, Mdx + 1, 1, Orb);

    // Record Coeficients Data
    // -----------------------

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[5]);
    if (timeinfo == 'r' || timeinfo == 'R')
    { strcat(fname_out, "_coef_realtime.dat");  }
    else
    { strcat(fname_out, "_coef_imagtime.dat"); }
    carr_txt(fname_out, mc->nc, C);

    // Record Energy Data
    // ------------------

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[5]);
    if (timeinfo == 'r' || timeinfo == 'R')
    { strcat(fname_out, "_E_realtime.dat");  }
    else
    { strcat(fname_out, "_E_imagtime.dat"); }
    carr_txt(fname_out, N, E);
    
    // Record Parameters Used
    // ----------------------

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[5]);
    if (timeinfo == 'r' || timeinfo == 'R')
    { strcat(fname_out, "_config_realtime.dat"); }
    else
    { strcat(fname_out, "_config_imagtime.dat"); }

    out_data = fopen(fname_out, "w");

    if (out_data == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_out); return -1; }

    fprintf(out_data, "%d %d %d %.15lf %.15lf %.10lf %d",
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
