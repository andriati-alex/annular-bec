#include <string.h>
#include "../include/MCTDHB_integrator.h"





/* ==========================================================================



   *  OBTAIN REAL OR IMAGINARY TIME EVOLUTION SOLUTION OF MCTDHB EQUATIONS  *



   REQUIRED FILES
   --------------------------------------------------------------------------

   (1)  setup/MC_fileId_orb.dat

        Text file with a matrix where the k-th column represent the k-th
        orbital. Thus the rows represent the values of these orbitals in
        discretized positions

   (2)  setup/MC_fileId_eq.dat

        A text file within the values of equation coefficients in the
        following order, separeted by spaces:

        (2-1) second order derivative
        (2-2) imag part of first order derivative(have no real part)
        (2-3) interaction strength (g)
        (2-5) Boundary conditions - Boolean

   (3)  setup/MC_fileId_config.dat

        A text file with position domain information over  which  the
        fileId_init was generated. The numbers are in following order

        (3-1) # of particles
        (3-2) # of orbitals
        (3-3) x_i
        (3-4) x_f
        (3-5) M the number of slices of size (xf - xi) / M





   COMMAND LINE ARGUMENTS
   --------------------------------------------------------------------------

   real/imag dt N fileId output_name method(optional) Nstates(optional)

        real/imag   -> real/Real/imag/Imag (the time)
        dt          -> the time step
        N           -> number of time steps to propagate
        fileId      -> Name to look up for files
        output_name -> name of file to write results
        method      -> Method to solve (default 1)
        Nstates     -> Find Ground state for multiple parameters (default 1)





   OBSERVATIONS
   --------------------------------------------------------------------------
   PS: method, if given, must be :

            1 Finite differences Crank-Nicolson with Sherman-Morrison
            2 Finite differences Crank-Nicolson with LU-decomposition
            3 Fast-Fourier Transform to compute derivatives

   Nstates is the number of lines in _config.dat and _eq.dat files
   to run multiple imaginary propagation  for  different  equation
   parameters set. Even if the file has several lines and  Nstates
   is not defined, the program will run once for the fisrt line.





   CALL
   --------------------------------------------------------------------------

   ./MCTDHB_time real/imag dt N fileId output_name method(optional)





   OUTPUT FILES
   --------------------------------------------------------------------------





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

    if ( tt / 60    > 0 )
    { mins  = tt / 60;    tt = tt % 60;    }

    printf("%d day(s) %d hour(s) %d minute(s)", days, hours, mins);
}





int main(int argc, char * argv[])
{

    omp_set_num_threads(omp_get_max_threads() / 2);
    mkl_set_num_threads(omp_get_max_threads() / 2);


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

    int
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
        xf,    // Domain of orbitals [xi, xf]
        dt,    // time step (both for real and imaginary)  
        real,  // real part of read data from file
        imag,  // imag part of read data from file
        a2,    // Term multiplying d2 / dx2
        inter, // contact interaction strength
        check, // to check norm
        * V;   // Potential computed in discretized positions

    double complex
        a1,        // Term multiplying d / dx
        checkDiag; // To check norm

    char
        timeinfo,
        fname_in[120],
        fname_out[120];
    
    FILE // pointer to file opened
        * eq_setup_file,
        * out_data;

    Carray
        dCdt,
        rho2,   // test time to setup 2-body density matris
        Ctest,  // estimate time based in on tiime-step
        C,      // Coeficients of superposition of Fock states
        to_int, // auxiliar to compute integration
        vir,    // virial at each time step (should be zero)
        E;      // Energy at each time step evolved

    Cmatrix
        rho,    // test time to setup 1-body density matrix
        Orbtest,
        Orb;    // Orb[k][j] give the value of k-th orbital at position j

    MCTDHBsetup
        mc;










    /* ==================================================================== *
     *                                                                      *
     *                    READ ALL DATA FROM SETUP FILES                    *
     *                                                                      *
     * ==================================================================== */



    printf("\n\n\n");
    printf("=========================================================\n\n");
    printf("Use MC_%s setup files to configure integrator\n\n", argv[4]);
    printf("=========================================================\n\n");



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
    { printf(" .... Found !\n"); }

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
    Orbtest = cmatDef(Morb, Mdx + 1);

    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[4]);
    strcat(fname_in, "_orb.dat");

    printf("\nLooking for %s ", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("\nERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ...... Found !\n"); }

    for (k = 0; k < Mdx + 1; k++)
    {
        for (s = 0; s < Morb; s++)
        {
            l = fscanf(eq_setup_file, " (%lf%lfj) ", &real, &imag);
            Orb[s][k] = real + I * imag;
            Orbtest[s][k] = real + I * imag;
        }
    }

    fclose(eq_setup_file); // finish the reading of file



    /* Setup Coeficients *
     * ----------------- */



    C = carrDef(NC(Npar, Morb));
    dCdt = carrDef(NC(Npar, Morb));
    Ctest = carrDef(NC(Npar, Morb));
    
    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[4]);
    strcat(fname_in, "_coef.dat");

    printf("\nLooking for %s ", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("\nERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ..... Found !\n"); }

    for (k = 0; k < NC(Npar, Morb); k++)
    {
        l = fscanf(eq_setup_file, " (%lf%lfj)", &real, &imag);
        C[k] = real + I * imag;
        Ctest[k] = real + I * imag;
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
    { printf(" ........ Found !\n"); }

    l = fscanf(eq_setup_file, "%lf %lf %lf %d",
                   &a2, &imag, &inter, &cyclic);

    fclose(eq_setup_file);

    a1 = 0 + imag * I;



    /* Setup Trap potential in discretized positions *
     * --------------------------------------------- */



    V = rarrDef(Mdx + 1);

    strcpy(fname_in, "setup/MC_");
    strcat(fname_in, argv[4]);
    strcat(fname_in, "_trap.dat");

    printf("\nLooking for %s ", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("\nERROR: impossible to open file %s\n", fname_in); return -1; }
    else
    { printf(" ..... Found !\n"); }

    for (k = 0; k < Mdx + 1; k++)
    {
        l = fscanf(eq_setup_file, " %lf", &real);
        V[k] = real;
    }

    fclose(eq_setup_file); // finish the reading of file



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



    if (argc > 6) { sscanf(argv[6], "%d", &method); }
    else          { method = 1;                     }

    if (method != 1 && method != 2 && method != 3)
    {
        printf("\n\n\tInvalid method Id ! Valid: 1 or 2 or 3\n\n");
        return -1;
    }

    // Print what it is going to do
    printf("\n\nMethod chosen : %d - ", method);
    switch (method)
    {
        case 1:
            printf("Crank-Nicolson SM and RK4 \n");
            break;
        case 2:
            printf("Crank-Nicolson LU and RK4 \n");
            break;
        case 3:
            printf("FFT and RK4 \n");
            break;
    }










    /* ==================================================================== *
     *                                                                      *
     *                      CHECK ORTHOGONAL CONDITION                      *
     *                                                                      *
     * ==================================================================== */



    printf("\n\n\n");
    printf("=========================================================\n\n");
    printf("Configuration done. Checking orthonormality\n\n");
    printf("=========================================================\n\n");


    mc = AllocMCTDHBdata(Npar, Morb, Mdx + 1, xi, xf, a2, inter, V, a1);

    to_int = carrDef(Mdx + 1);

    E = carrDef(N + 1);   // to store energy
    vir = carrDef(N + 1); // check consistency by Virial Theorem



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

    if (abs(creal(checkDiag) - Morb) > 1E-8 || cimag(checkDiag) > 1E-8)
    {
        printf("\n\n\tOrbitals are not normalized to 1 !\n"); return -1;
    }


    /* Check normalization of coeficients *
     * ---------------------------------- */


    if ( abs(carrMod2(NC(Npar, Morb), C) - 1) > 1E-9 )
    {
        printf("\n\n\tCoeficients norm is not 1 !\n");
        return -1;
    }





    /* ==================================================================== *
     *                                                                      *
     *                           SOME TIME REQUIRED                         *
     *                                                                      *
     * ==================================================================== */

    rho2 = carrDef(Morb * Morb * Morb * Morb);
    rho  = cmatDef(Morb, Morb);

    start = omp_get_wtime();
    OBrho(Npar, Morb, mc->NCmat, mc->IF, C, rho);
    time_used = (double) (omp_get_wtime() - start);

    printf("\n\ntime to setup rho = %.3lf", time_used);
    
    start = omp_get_wtime();
    TBrho(Npar, Morb, mc->NCmat, mc->IF, C, rho2);
    time_used = (double) (omp_get_wtime() - start);

    printf("\n\ntime to setup rho2 = %.3lf", time_used);

    start = omp_get_wtime();
    SetupHo(Morb, Mdx + 1, Orb, mc->dx, a2, a1, V, rho);
    time_used = (double) (omp_get_wtime() - start);
    
    printf("\n\ntime to setup Ho = %.3lf", time_used);

    start = omp_get_wtime();
    SetupHint(Morb, Mdx + 1, Orb, mc->dx, inter, rho2);
    time_used = (double) (omp_get_wtime() - start);

    printf("\n\ntime to setup Hint = %.3lf", time_used);

    start = omp_get_wtime();
    applyHconf(mc->Npar, mc->Morb, mc->NCmat, mc->IF, C, rho, rho2, dCdt);
    time_used = (double) (omp_get_wtime() - start);

    printf("\n\ntime to apply Many-Body H = %.3lf\n\n", time_used);










    /* ==================================================================== *
     *                                                                      *
     *       DIAGONALIZE HAMILTONIAN IN THE GIVEN BASIS BEFORE START        *
     *                                                                      *
     * ==================================================================== */





    if ( dt > 5 * dx * dx && method < 3)
    {
        printf("\n\nWARNING : step too large to maintain stability");
        printf(" in finite-differences methods.\n\n");
    }



    if ( NC(Npar, Morb) < 400 )
    {
        E[0] = LanczosGround( NC(Npar,Morb)/2, mc, Orb, C );
        // Renormalize coeficients
        renormalizeVector( NC(Npar,Morb), C, 1.0);
    } else
    {
        E[0] = LanczosGround( 200, mc, Orb, C );
        // Renormalize coeficients
        renormalizeVector( NC(Npar,Morb), C, 1.0);
    }



    // Test if time step is good
    if ( dt > 0.2 / creal(E[0]) )
    {
        printf("\n\n\n\t!   Too big time step   !");
        printf("\n\nTry something < %.10lf\n\n", 0.15 / creal(E[0]));
        exit(EXIT_FAILURE);
    } else
    {
        if ( dt < 0.001 / creal(E[0]) )
        {
            printf("\n\n\n\t!   Too small time step   !");
            printf("\n\nTry something > %.10lf\n\n", 0.003 / creal(E[0]));
            exit(EXIT_FAILURE);
        }
    }



    // Test if final time is good
    if ( N * dt < 20 / creal(E[0]) )
    {
        printf("\n\n\n\t!   Final step(%.3lf) too small   !", N * dt);
        printf("\n\nTry N * dt > %.3lf\n\n", 21 / creal(E[0]));
        exit(EXIT_FAILURE);
    }





    /* ==================================================================== *
     *                                                                      *
     *                          CALL THE INTEGRATOR                         *
     *                                                                      *
     * ==================================================================== */



    printf("\n\n\n");
    printf("=========================================================\n\n");
    printf("Start Integration in pure imaginary time\n\n");
    printf("=========================================================\n\n");

    // setup filename to store solution
    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[5]);

    switch (method)
    {

        case 1:

            start = omp_get_wtime();
            MC_IMAG_RK4_CNSMRK4(mc, Orbtest, Ctest, E, vir, dt, 1, cyclic);
            time_used = (double) (omp_get_wtime() - start);

            printf("\n\nTime to do 1 step: %.3lf seconds", time_used);
            printf("\nTotal time estimated: ");
            TimePrint(time_used * N);

            free(Ctest);
            cmatFree(Morb, Orbtest);

            printf("\n\n");

            // Start Evolution
            start = omp_get_wtime();
            MC_IMAG_RK4_CNSMRK4(mc, Orb, C, E, vir, dt, N, cyclic);
            time_used = (double) (omp_get_wtime() - start);

            break;

        case 2:

            start = omp_get_wtime();
            MC_IMAG_RK4_CNLURK4(mc, Orbtest, Ctest, E, vir, dt, 1, cyclic);
            time_used = (double) (omp_get_wtime() - start);

            printf("\n\nTime to do 1 step: %.3lf seconds", time_used);
            printf("\nTotal time estimated: ");
            TimePrint(time_used * N);

            free(Ctest);
            cmatFree(Morb, Orbtest);

            printf("\n\n");

            // Start Evolution
            start = omp_get_wtime();
            MC_IMAG_RK4_CNLURK4(mc, Orb, C, E, vir, dt, N, cyclic);
            time_used = (double) (omp_get_wtime() - start);

            break;

        case 3:

            start = omp_get_wtime();
            MC_IMAG_RK4_FFTRK4(mc, Orbtest, Ctest, E, vir, dt, 1);
            time_used = (double) (omp_get_wtime() - start);

            printf("\n\nTime to do 1 step: %.3lf seconds", time_used);
            printf("\nTotal time estimated: ");
            TimePrint(time_used * N);

            free(Ctest);
            cmatFree(Morb, Orbtest);

            printf("\n\n");

            // Start Evolution
            start = omp_get_wtime();
            MC_IMAG_RK4_FFTRK4(mc, Orb, C, E, vir, dt, N);
            time_used = (double) (omp_get_wtime() - start);

            break;
    }










    /* ==================================================================== *
     *                                                                      *
     *                              Record Data                             *
     *                                                                      *
     * ==================================================================== */


    // Record data in case of imaginary time
    // ---------------------------------------------

    if (timeinfo == 'i' || timeinfo == 'I')
    {
        // record orbital data

        strcat(fname_out, "_orb_imagtime.dat");
        cmat_txt(fname_out, Morb, 1, Mdx + 1, 1, Orb);

        // Record Coeficients Data

        strcpy(fname_out, "../mctdhb_data/");
        strcat(fname_out, argv[5]);
        strcat(fname_out, "_coef_imagtime.dat");

        carr_txt(fname_out, mc->nc, C);

        // Record Energy Data in case of imaginary time

        strcpy(fname_out, "../mctdhb_data/");
        strcat(fname_out, argv[5]);
        strcat(fname_out, "_E_imagtime.dat");

        carr_txt(fname_out, N + 1, E);
        
        strcpy(fname_out, "../mctdhb_data/");
        strcat(fname_out, argv[5]);
        strcat(fname_out, "_virial_imagtime.dat");

        carr_txt(fname_out, N + 1, vir);
    }



    // Record Parameters Used

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[5]);

    if (timeinfo == 'r' || timeinfo == 'R')
    { strcat(fname_out, "_conf_realtime.dat"); }
    else
    { strcat(fname_out, "_conf_imagtime.dat"); }

    out_data = fopen(fname_out, "w");

    if (out_data == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n", fname_out);
        return -1;
    }

    fprintf(out_data, "%d %d %d %.15lf %.15lf %.15lf %.15lf %.15lf",
            Npar, Morb, Mdx, xi, xf, mc->a2, cimag(mc->a1), mc->inter);

    fclose(out_data);
    
    
    
    // Record Trap potential

    strcpy(fname_out, "../mctdhb_data/");
    strcat(fname_out, argv[5]);

    if (timeinfo == 'r' || timeinfo == 'R')
    { strcat(fname_out, "_trap_realtime.dat"); }
    else
    { strcat(fname_out, "_trap_imagtime.dat"); }

    rarr_txt(fname_out, Mdx + 1, mc->V);



    /* ==================================================================== *
     *                                                                      *
     *                      FINISH UP - RELEASE MEMORY                      *
     *                                                                      *
     * ==================================================================== */



    EraseMCTDHBdata(mc);
    free(C);
    free(E);
    free(vir);
    free(dCdt);
    free(rho2);
    free(to_int);
    cmatFree(Morb, Orb);
    cmatFree(Morb, rho);

    printf("\n\n");
    return 0;
}
