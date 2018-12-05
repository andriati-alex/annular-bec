#include <string.h>
#include "../include/MCTDHB_integrator.h"
#include "../include/linear_potential.h"





/* ==========================================================================



   *  OBTAIN REAL OR IMAGINARY TIME EVOLUTION SOLUTION OF MCTDHB EQUATIONS  *



   REQUIRED FILES
   --------------------------------------------------------------------------

   (1)  setup/MC_fileId_orb.dat

        Text file with a matrix where the k-th column represent the k-th
        orbital. Thus the rows represent the values of these orbitals in
        discretized positions and each column is then an orbital.

   (2)  setup/MC_fileId_eq.dat

        A text file within the values of equation coefficients organized
        by columns, separeted by spaces:

        Col 1 - second order derivative
        Col 2 - imag part of first order derivative(have no real part)
        Col 3 - interaction strength (g)

        Each line is used as one possible setup to find a ground state

   (3)  setup/MC_fileId_config.dat

        A text file with position/time domain information over  which
        the orbitals and coefficients were generated. The numbers are
        in columns separeted by spaces, in each line as follows:

        Col 1 - # of particles (may vary along lines)
        Col 2 - # of orbitals
        Col 5 - (Mpos) # of slices of size (xf - xi) / Mpos
        Col 3 - x_i
        Col 4 - x_f
        Col 5 - dt (may vary along lines)
        Col 6 - # of time-steps (may vary along lines)





   COMMAND LINE ARGUMENTS
   --------------------------------------------------------------------------

   time_domain cyclic fileId_in FileId_out method(optional) Nstates(optional)

        time_domain -> real/Real/imag/Imag
        cyclic      -> 1 or 0 (boolean)
        fileId_in   -> Name to look up for files
        fileId_out  -> name of file to write results
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
   is not defined, the program will run once for the fisrt line. A
   reasonable step for the varying parameters must be chosen since
   the program take the initial orbitals those from  the  previous
   propagation. If some parameter vary abruptely  then this choice
   may not be the ideal and problems in convergence are likely  to
   occur.

   NOTE THAT  # OF ORBITALS AND POSITION DOMAIN DISCRETIZATION ARE
   NOT ALLOWED TO CHANGE BETWEEN LINES.

   Every line the # of particles change, then the file '_coef.dat'
   is the opened again to read the coefficients. The program then
   assume the new vector of coefficients is concatenated to   the
   previous one. For example, if the # of particles change  three
   times among the lines in _conf.dat file, lets say N1 N2 and N3
   thus is read from the file  NC( N1 , Morb ) elements taken  as
   initial condition, and when it changes to N2 the program  read
   more NC( N2 , Morb ) elements from the file taken again as the
   new initial condition and so on.

   If the # of particles does not vary between two lines then the
   new initial conditions is taken as the last one, just like  is
   done for orbitals explained above





   CALL
   --------------------------------------------------------------------------

   ./MCTDHB_time time_domain cyclic fileId_in fileId_out
   ./MCTDHB_time time_domain cyclic fileId_in fileId_out method
   ./MCTDHB_time time_domain cyclic fileId_in fileId_out method Nstates





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










void SaveConf(char timeinfo, char fnameIn [], char fnameOut [], int Nlines)
{

    int
        i,
        N,
        k,
        Npar,
        Morb,
        Mdx;

    double
        p1,
        p2,
        p3,
        xi,
        xf,
        dt,
        a2,
        g,
        imag;

    char
        Vname[30],
        fname [120];

    FILE
        * confFileIn,
        * eqFileIn,
        * confFileOut;

    strcpy(fname, "setup/MC_");
    strcat(fname, fnameIn);
    strcat(fname, "_conf.dat");

    confFileIn = fopen(fname, "r");
    if (confFileIn == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        return;
    }

    strcpy(fname, "setup/MC_");
    strcat(fname, fnameIn);
    strcat(fname, "_eq.dat");

    eqFileIn = fopen(fname, "r");
    if (eqFileIn == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        return;
    }
    
    strcpy(fname, "../mctdhb_data/");
    strcat(fname, fnameOut);

    if (timeinfo == 'r' || timeinfo == 'R')
    { strcat(fname, "_conf_realtime.dat"); }
    else
    { strcat(fname, "_conf_imagtime.dat"); }

    confFileOut = fopen(fname, "w");
    if (confFileOut == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        return;
    }
    
    // Read data and write in out file (transfer)

    for (i = 0; i < Nlines; i++)
    {
        k = fscanf(eqFileIn, "%lf %lf %lf %s %lf %lf %lf",
                   &a2, &imag, &g, Vname, &p1, &p2, &p3);

        k = fscanf(confFileIn, "%d %d %d %lf %lf %lf %d",
                   &Npar, &Morb, &Mdx, &xi, &xf, &dt, &N);

        fprintf(confFileOut, "%d %d %d %.15lf %.15lf %.15lf %d ",
                Npar, Morb, Mdx, xi, xf, dt, N);

        fprintf(confFileOut, "%.15lf %.15lf %.15lf ", a2, imag, g);

        fprintf(confFileOut, "%s %.15lf %.15lf %.15lf", Vname, p1, p2, p3);

        fprintf(confFileOut, "\n");
    }

    fclose(eqFileIn); fclose(confFileIn); fclose(confFileOut);

}










MCTDHBsetup SetupData(FILE * paramFile, FILE * confFile, Rarray x,
            double * dt, int * N)
{

    int
        k,
        Mdx,
        Npar,
        Morb;

    char
        Vname[30];

    double
        p1,
        p2,
        p3,
        xi,
        xf,
        dx,
        a2,
        imag,
        inter;

    double complex
        a1;

    Rarray
        V;



    // Setup spatial, time, num of particles and orbitals
    // --------------------------------------------------

    k = fscanf(confFile, "%d %d %d %lf %lf %lf %d",
               &Npar, &Morb, &Mdx, &xi, &xf, dt, N);

    dx = (xf - xi) / Mdx;

    V = rarrDef(Mdx + 1);

    printf("\n\nConfiguration Done\n");

    printf("\n\t# of Particles: %d", Npar);
    printf("\n\t# of Orbitals: %d", Morb);
    printf("\n\t# of possible configurations: %d\n", NC(Npar, Morb));



    // Setup Equation parameters
    // -------------------------

    k = fscanf(paramFile, "%lf %lf %lf %s %lf %lf %lf",
               &a2, &imag, &inter, Vname, &p1, &p2, &p3);

    GetPotential(Mdx + 1, Vname, x, V, p1, p2, p3);

    a1 = 0 + imag * I;

    return AllocMCTDHBdata(Npar, Morb, Mdx + 1, xi, xf, a2, inter, V, a1);

}










int main(int argc, char * argv[])
{

    omp_set_num_threads(omp_get_max_threads() / 2);
    mkl_set_num_threads(omp_get_max_threads() / 2);


    if (argc < 5 || argc > 7)
    {
        printf("\nInvalid Number of command line arguments, ");
        printf("expected: 4,5 or 6.\n\n");
        return -1;
    }

    /* ==================================================================== *
     *                                                                      *
     *                      VARIABLES AND WHAT THEY DO                      *
     *                                                                      *
     * ==================================================================== */

    int
        i,
        k,
        l,
        s;

    int
        N,    // # of time steps to evolve the system
        Mdx,  // # of divisions in space (# of points - 1)
        Npar, // # of particles
        Morb, // # of orbitals
        cyclic,
        Nlines,
        method;

    double
        start, // trigger to measure time
        end,   // trigger to finish time per execution
        time_used, // total time used
        dx,
        xi,
        xf,    // Domain of orbitals [xi, xf]
        dt,    // time step (both for real and imaginary)  
        real,  // real part of read data from file
        imag,  // imag part of read data from file
        check, // to check norm
        * x;   // discretized positions

    double complex
        checkDiag; // To check norm

    char
        timeinfo,
        strnum[30],
        fname[120];
    
    FILE // pointer to file opened
        * E_file,
        * coef_file,
        * orb_file,
        * confFile,
        * paramFile;

    Carray
        C,      // Coeficients of superposition of Fock states
        to_int, // auxiliar to compute integration
        vir,    // virial at each time step (should be zero)
        E;      // Energy at each time step evolved

    Cmatrix
        Orb;    // Orb[k][j] give the value of k-th orbital at position j

    MCTDHBsetup
        mc;










    /* ====================================================================
                          CONFIGURE TYPE OF INTEGRATION
       ==================================================================== */

    timeinfo = argv[1][0]; // character i for imaginary or r for real time

    if (timeinfo != 'r' && timeinfo != 'R')
    {
        if (timeinfo != 'i' && timeinfo != 'I')
        {
            printf("\n\n\t!   Invalid first argument   !\n");
            return -1;
        }
    }

    sscanf(argv[2], "%d", &cyclic);

    if (argc > 5) { sscanf(argv[5], "%d", &method); }
    else          { method = 1;                     }

    if (method != 1 && method != 2 && method != 3)
    {
        printf("\n\n\tInvalid method Id ! Valid: 1 or 2 or 3\n\n");
        return -1;
    }

    if (argc == 7) { sscanf(argv[6], "%d", &Nlines); }
    else           { Nlines = 1;                     }

    // Print what it is going to do
    printf("\n\nMethod chosen : %d - ", method);
    switch (method)
    {
        case 1:
            printf("Crank-Nicolson SM and RK4");
            break;
        case 2:
            printf("Crank-Nicolson LU and RK4");
            break;
        case 3:
            printf("FFT and RK4");
            break;
    }










    /* ====================================================================
                   CONFIGURE TRAP POTENTIAL AND DISCRETIZATION
       ==================================================================== */

    printf("\n\n\n");
    printf("=========================================================\n\n");
    printf("Using MC_%s setup files to configure integrator\n\n", argv[3]);
    printf("=========================================================\n\n");



    // Read number of discrete positions in spatial domain
    // ---------------------------------------------------

    strcpy(fname, "setup/MC_");
    strcat(fname, argv[3]);
    strcat(fname, "_conf.dat");

    printf("\n\nLooking for %s", fname);

    confFile = fopen(fname, "r");

    if (confFile == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        return -1;
    } else
    {
        printf(" .... Found !\n");
    }

    // Get the third value - number of domain slices
    k = fscanf(confFile, "%d ", &Mdx);
    k = fscanf(confFile, "%d ", &Mdx);
    k = fscanf(confFile, "%d ", &Mdx);
    // Get the domain limit
    k = fscanf(confFile, "%lf ", &xi);
    k = fscanf(confFile, "%lf ", &xf);

    dx = (xf - xi) / Mdx;
    
    x  = rarrDef(Mdx + 1);
    
    rarrFillInc(Mdx + 1, xi, dx, x);

    fclose(confFile);



    // Record the setup file to further analisys
    SaveConf(timeinfo, argv[3], argv[4], Nlines);









    
    /* ====================================================================
                         OPEN FILES TO SETUP THE PROBLEM
       ==================================================================== */

    strcpy(fname, "setup/MC_");
    strcat(fname, argv[3]);
    strcat(fname, "_conf.dat");

    confFile = fopen(fname, "r");



    strcpy(fname, "setup/MC_");
    strcat(fname, argv[3]);
    strcat(fname, "_eq.dat");

    printf("\nLooking for %s", fname);

    paramFile = fopen(fname, "r");

    if (paramFile == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        return -1;
    } else
    {
        printf(" ...... Found !\n");
    }



    strcpy(fname, "setup/MC_");
    strcat(fname, argv[3]);
    strcat(fname, "_orb.dat");

    printf("\nLooking for %s ", fname);

    orb_file = fopen(fname, "r");

    if (orb_file == NULL)  // impossible to open file
    {
        printf("\n\nERROR: impossible to open file %s\n", fname);
        return -1;
    } else
    {
        printf(" .... Found !\n");
    }



    strcpy(fname, "setup/MC_");
    strcat(fname, argv[3]);
    strcat(fname, "_coef.dat");

    printf("\nLooking for %s ", fname);

    coef_file = fopen(fname, "r");

    if (coef_file == NULL)  // impossible to open file
    {
        printf("\nERROR: impossible to open file %s\n", fname);
        return -1;
    } else
    {
        printf(" ... Found !\n");
    }



    // open file to write energy values
    // ----------------------------------------------------------------
    strcpy(fname, "../mctdhb_data/");
    strcat(fname, argv[4]);
    strcat(fname, "_energy_imagtime.dat");

    E_file = fopen(fname, "w");

    if (E_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }










    /* ====================================================================
           READ DATA TO SETUP EQUATION PARAMETERS AND INITIAL CONDITIONS
       ==================================================================== */

    mc = SetupData(paramFile, confFile, x, &dt, &N);

    Morb = mc->Morb;
    Npar = mc->Npar;



    // Setup orbitals
    // --------------

    Orb = cmatDef(Morb, Mdx + 1);

    for (k = 0; k < Mdx + 1; k++)
    {
        for (s = 0; s < Morb; s++)
        {
            l = fscanf(orb_file, " (%lf%lfj) ", &real, &imag);
            Orb[s][k] = real + I * imag;
        }
    }

    // orbitals are not read again
    fclose(orb_file);



    // Setup Coeficients
    // -----------------

    C = carrDef(NC(Npar, Morb));

    for (k = 0; k < NC(Npar, Morb); k++)
    {
        l = fscanf(coef_file, " (%lf%lfj)", &real, &imag);
        C[k] = real + I * imag;
    }










    /* ====================================================================
                       CHECK ORTHOGONALITY AND NORMALIZATION
       ==================================================================== */

    printf("\n\n\n");
    printf("=========================================================\n\n");
    printf("Configuration done. Checking orthonormality\n\n");
    printf("=========================================================\n\n");


    to_int = carrDef(Mdx + 1);



    // Check if off-diagonal elements are zero
    // ---------------------------------------

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
            check = check + cabs(Csimps(Mdx + 1, to_int, mc->dx));
        }
    }

    if (check > 1E-8)
    {
        printf("\n\n\t!   ORBITALS ARE NOT ORTHOGONAL   !\n");
        return -1;
    }



    // Check if Diagonal elements sum up to Morb
    // -----------------------------------------

    checkDiag = 0;
    for (k = 0; k < Morb; k++)
    {
        for (s = 0; s < Mdx + 1; s++)
        {
            to_int[s] = conj(Orb[k][s]) * Orb[k][s];
        }
        checkDiag = checkDiag + Csimps(Mdx + 1, to_int, mc->dx);
    }

    if (abs(creal(checkDiag) - Morb) > 1E-8 || cimag(checkDiag) > 1E-8)
    {
        printf("\n\n\t!   ORBITALS DO NOT HAVE NORM = 1   !\n");
        return -1;
    }



    // Check normalization of coeficients
    // ----------------------------------

    if ( abs(carrMod2(NC(Npar, Morb), C) - 1) > 1E-9 )
    {
        printf("\n\n\t!   COEFFICIENTS DO NOT HAVE NORM = 1   !\n");
        return -1;
    }










    /* ====================================================================
             DIAGONALIZE HAMILTONIAN IN THE GIVEN BASIS BEFORE START
       ==================================================================== */

    E   = carrDef(N + 1); // to store energy
    vir = carrDef(N + 1); // check consistency by Virial Theorem

    if ( dt > 5 * dx * dx && method < 3)
    {
        printf("\n\nWARNING : time-step too large to maintain stability");
        printf(" in finite-differences methods.\n\n");
    }



    // Estimate energy by diagonalization Restrict the number of iterations
    // in lanczos routine to avoid massive memory usage. Try to use 200
    // iterations unless either it exceeds half of the dimension of
    // configuration space or if it exceeds a memory Threshold.

    if (200 * NC(Npar, Morb) < 5E7)
    {
        if (NC(Npar, Morb) / 2 < 200) k = NC(Npar, Morb) / 2;
        else                          k = 200;
    }
    else k = 5E7 / NC(Npar, Morb);

    E[0] = LanczosGround( k, mc, Orb, C );
    // Renormalize coeficients
    renormalizeVector(NC(Npar, Morb), C, 1.0);

    printf("\n=====================================================");
    printf("==========================\n\n");
    printf("\tDiagonalization Done E = %.5E", creal(E[0]));
    printf("\n\t# of lanczos iterations used = %d", k);
    printf("\n\n===================================================");
    printf("============================");


    // Test if time step is good
    if ( dt > 0.1 / creal(E[0]) )
    {
        printf("\n\n\n\t!   WARNING : Large time step   !");
        printf("\n\nTry something < %.10lf\n\n", 0.09 / creal(E[0]));
    } else
    {
        if ( dt < 0.001 / creal(E[0]) )
        {
            printf("\n\n\n\t!   WARNING : Too small time step   !");
            printf("\n\nTry something > %.10lf\n\n", 0.003 / creal(E[0]));
        }
    }



    // Test if final time is good
    if ( N * dt < 15 / creal(E[0]) )
    {
        printf("\n\n\n\t!   WARNING : Final step(%.3lf) too small   !", N*dt);
        printf("\n\nTry N * dt > %.3lf\n\n", 17 / creal(E[0]));
    }










    /* ====================================================================
                                CALL THE INTEGRATOR 
     * ==================================================================== */

    printf("\n\n\n");
    printf("=========================================================\n\n");
    printf("Start Integration #%d\n\n", 1);
    printf("=========================================================\n\n");

    // setup filename to record solution
    strcpy(fname, "../mctdhb_data/");
    strcat(fname, argv[4]);
    strcat(fname, "_line-1");

    switch (method)
    {

        case 1:

            start = omp_get_wtime();
            s = MC_IMAG_RK4_CNSMRK4(mc, Orb, C, E, vir, dt, N, cyclic);
            time_used = (double) (omp_get_wtime() - start);
            printf("\n\nTime taken in integration #%d : %lf", 1, time_used);
            printf(" = "); TimePrint(time_used);

            break;

        case 2:

            start = omp_get_wtime();
            s = MC_IMAG_RK4_CNLURK4(mc, Orb, C, E, vir, dt, N, cyclic);
            time_used = (double) (omp_get_wtime() - start);
            printf("\n\nTime taken in integration #%d : %lf", 1, time_used);
            printf(" = "); TimePrint(time_used);

            break;

        case 3:

            start = omp_get_wtime();
            s = MC_IMAG_RK4_FFTRK4(mc, Orb, C, E, vir, dt, N);
            time_used = (double) (omp_get_wtime() - start);
            printf("\n\nTime taken in integration #%d : %lf", 1, time_used);
            printf(" = "); TimePrint(time_used);

            break;
    }



    // Record data
    // ---------------------------------------------

    if (timeinfo == 'i' || timeinfo == 'I')
    {
        // record orbital data

        strcat(fname, "_orb_imagtime.dat");
        cmat_txt(fname, Morb, 1, Mdx + 1, 1, Orb);

        // Record Coeficients Data

        strcpy(fname, "../mctdhb_data/");
        strcat(fname, argv[4]);
        strcat(fname, "_line-1");
        strcat(fname, "_coef_imagtime.dat");

        carr_txt(fname, mc->nc, C);

        // Record Energy
        /*
        strcpy(fname, "../mctdhb_data/");
        strcat(fname, argv[4]);
        strcat(fname, "_line-1");
        strcat(fname, "_E_imagtime.dat");

        carr_txt(fname, s, E);

        strcpy(fname, "../mctdhb_data/");
        strcat(fname, argv[4]);
        strcat(fname, "_line-1");
        strcat(fname, "_virial_imagtime.dat");

        carr_txt(fname, s, vir);
        */

        fprintf(E_file, "%.10E\n", creal(E[s-1]));
    }











/** If either the _conf.dat or _eq.dat file have more  than  one  line
    and the Nlines parameter is greater than 1 it read the next config
    from files. **/

    for (i = 1; i < Nlines; i++)
    {

        // number of line reading in _conf.dat and _eq.dat files
        sprintf(strnum, "%d", i + 1);

        // release old data
        EraseMCTDHBdata(mc);

        // setup new parameters
        mc = SetupData(paramFile, confFile, x, &dt, &N);

        if (Npar != mc->Npar)
        {

        //  Setup Coeficients from file because the number
        //  of particles has changed

            free(C);

            Npar = mc->Npar;

            C = carrDef(NC(Npar, Morb));

            for (k = 0; k < NC(Npar, Morb); k++)
            {
                l = fscanf(coef_file, " (%lf%lfj)", &real, &imag);
                C[k] = real + I * imag;
            }
        }

        // Use orbitals from previous result



        /* ================================================================
                         CHECK ORTHOGONALITY AND NORMALIZATION
           ================================================================ */

        printf("\n\n\n");
        printf("=========================================================\n\n");
        printf("Configuration done. Checking orthonormality\n\n");
        printf("=========================================================\n\n");

        // Check if off-diagonal elements are zero

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
                check = check + cabs(Csimps(Mdx + 1, to_int, mc->dx));
            }
        }

        if (check > 1E-9)
        {
            printf("\n\n\t!   ORBITALS ARE NOT ORTHOGONAL   !\n");
            return -1;
        }

        // Check if Diagonal elements sum up to Morb

        checkDiag = 0;
        for (k = 0; k < Morb; k++)
        {
            for (s = 0; s < Mdx + 1; s++)
            {
                to_int[s] = conj(Orb[k][s]) * Orb[k][s];
            }
            checkDiag = checkDiag + Csimps(Mdx + 1, to_int, mc->dx);
        }

        if (abs(creal(checkDiag) - Morb) > 1E-8 || cimag(checkDiag) > 1E-8)
        {
            printf("\n\n\t!   ORBITALS DO NOT HAVE NORM = 1   !\n");
            return -1;
        }

        // Check normalization of coeficients

        if ( abs(carrMod2(NC(Npar, Morb), C) - 1) > 1E-9 )
        {
            printf("\n\n\t!   COEFFICIENTS DO NOT HAVE NORM = 1   !\n");
            return -1;
        }



        /* =================================================================
                DIAGONALIZE HAMILTONIAN IN THE GIVEN BASIS BEFORE START
           ================================================================= */

        free(E);
        free(vir);

        E   = carrDef(N + 1); // to store energy
        vir = carrDef(N + 1); // check consistency by Virial Theorem

        if ( dt > 5 * dx * dx && method < 3)
        {
            printf("\n\nWARNING : time-step too large to maintain stability");
            printf(" in finite-differences methods.\n\n");
        }

        if (200 * NC(Npar, Morb) < 5E7)
        {
            if (NC(Npar, Morb) / 2 < 200) k = NC(Npar, Morb) / 2;
            else                          k = 200;
        }
        else k = 5E7 / NC(Npar, Morb);

        E[0] = LanczosGround( k, mc, Orb, C );
        // Renormalize coeficients
        renormalizeVector(NC(Npar, Morb), C, 1.0);

        printf("\n=====================================================");
        printf("==========================\n\n");
        printf("\tDiagonalization Done E = %.5E", creal(E[0]));
        printf("\n\t# of lanczos iterations used = %d", k);
        printf("\n\n===================================================");
        printf("============================");

        // Test if time step is good
        if ( dt > 0.1 / creal(E[0]) )
        {
            printf("\n\n\n\t!   WARNING : Large time step   !");
            printf("\n\nTry something < %.10lf\n\n", 0.09 / creal(E[0]));
        } else
        {
            if ( dt < 0.001 / creal(E[0]) )
            {
                printf("\n\n\n\t!   WARNING : Too small time step   !");
                printf("\n\nTry something > %.10lf\n\n", 0.003 / creal(E[0]));
            }
        }

        // Test if final time is good
        if ( N * dt < 15 / creal(E[0]) )
        {
            printf("\n\n\n\t! WARNING : Final step(%.3lf) too small !", N*dt);
            printf("\n\nTry N * dt > %.3lf\n\n", 17 / creal(E[0]));
        }



        /* ================================================================
                                 CALL THE INTEGRATOR 
         * ================================================================ */

        printf("\n\n\n");
        printf("=======================================================\n\n");
        printf("Start Integration #%d\n\n", i + 1);
        printf("=======================================================\n\n");

        // setup filename to store solution
        strcpy(fname, "../mctdhb_data/");
        strcat(fname, argv[4]);
        strcat(fname, "_line-");
        strcat(fname, strnum);

        switch (method)
        {

            case 1:

                start = omp_get_wtime();
                s = MC_IMAG_RK4_CNSMRK4(mc, Orb, C, E, vir, dt, N, cyclic);
                end = (double) (omp_get_wtime() - start);
                time_used += end;
                printf("\n\nTime taken in execution #%d : %.1lf", i+1, end);
                printf(" = "); TimePrint(end);

                break;

            case 2:

                start = omp_get_wtime();
                s = MC_IMAG_RK4_CNLURK4(mc, Orb, C, E, vir, dt, N, cyclic);
                end = (double) (omp_get_wtime() - start);
                time_used += end;
                printf("\n\nTime taken in execution #%d : %.1lf", i+1, end);
                printf(" = "); TimePrint(end);

                break;

            case 3:

                start = omp_get_wtime();
                s = MC_IMAG_RK4_FFTRK4(mc, Orb, C, E, vir, dt, N);
                end = (double) (omp_get_wtime() - start);
                time_used += end;
                printf("\n\nTime taken in execution #%d : %.1lf", i+1, end);
                printf(" = "); TimePrint(end);

                break;
        }



        // Record data
        // ---------------------------------------------

        if (timeinfo == 'i' || timeinfo == 'I')
        {
            // record orbital data

            strcat(fname, "_orb_imagtime.dat");
            cmat_txt(fname, Morb, 1, Mdx + 1, 1, Orb);

            // Record Coeficients Data

            strcpy(fname, "../mctdhb_data/");
            strcat(fname, argv[4]);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_coef_imagtime.dat");

            carr_txt(fname, mc->nc, C);

            // Record Energy
            /*
            strcpy(fname, "../mctdhb_data/");
            strcat(fname, argv[4]);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_E_imagtime.dat");

            carr_txt(fname, s, E);
            
            strcpy(fname, "../mctdhb_data/");
            strcat(fname, argv[4]);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_virial_imagtime.dat");

            carr_txt(fname, s, vir);
            */

            fprintf(E_file, "%.10E\n", creal(E[s-1]));
        }

    }





    /* ====================================================================
                                  RELEASE MEMORY
     * ==================================================================== */

    fclose(E_file);
    fclose(confFile);
    fclose(paramFile);
    fclose(coef_file);

    EraseMCTDHBdata(mc);
    free(x);
    free(C);
    free(E);
    free(vir);
    free(to_int);
    cmatFree(Morb, Orb);

    printf("\n\nTotal time taken: %.1lf = ", time_used);
    TimePrint(time_used);
    
    printf("\n\nAverage time per state: %.1lf = ", time_used / Nlines);
    TimePrint(time_used / Nlines);

    printf("\n\n");
    return 0;
}
