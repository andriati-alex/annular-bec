#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../include/GP_realtime_integrator.h"
#include "../include/GP_imagtime_integrator.h"
#include "../include/linear_potential.h"





/* =========================================================================
 *
 *
 * OBTAIN REAL OR IMAGINARY TIME EVOLUTION SOLUTION OF
 * GROSS-PITAEVSKII EQUATION WITH  DIRAC DELTA BARRIER
 *
 *
 *
 *
 *
 *
 * REQUIRED FILES
 * -------------------------------------------------------------------------
 *
 * setup/fileId_init.dat
 *
 *      is a text file with an array of complex numbers of size M + 1
 *      where M is the number of space slices of the domain. Must  be
 *      formatted according to numpy.savetxt function.
 *
 * setup/fileId_eq.dat
 *
 *      is a text file with values of equation coefficients in the
 *      following order:
 *
 *      (1) second order derivative
 *      (2) imag part of first order derivative(have no real part)
 *      (3) interaction strength
 *      (4) Delta Barrier Strength
 *      (5) Boundary conditions - Boolean
 *
 * setup/fileId_domain.dat
 *
 *      a text file with position domain information over  which  the
 *      fileId_init was generated. The numbers are in following order
 *
 *      (1) x_i
 *      (2) x_f
 *      (3) M the number of slices in domain of size (xf - xi) / M
 *
 *
 *
 *
 *
 * COMMAND LINE ARGUMENTS
 * -------------------------------------------------------------------------
 * 
 * dt N fileId method(optional)
 *
 *      dt     -> the time step
 *      N      -> number of time steps to evolve
 *      fileId -> the prefix of required file name
 *      method -> Method to solve the linear part(optional)
 *      timeId -> Real or imaginary time (real if not passed)
 *
 * PS: method, if given, must be:
 *
 *      1 Crank-Nicolson with Sherman-Morrison
 *      2 Crank-Nicolson with special LU-decomposition
 *      3 4th order Runge-Kutta for non linear part
 *      4 Fast-Fourier Transform
 *
 *
 *
 *
 *
 * CALL
 * -------------------------------------------------------------------------
 *
 * ./time_evolution dt N fileId method(optional)
 *
 *
 *
 *
 *
 * OUTPUT FILES
 * -------------------------------------------------------------------------
 *
 * ../gp_data/fileId_time.dat
 *
 *      a text file of matrix-like format with solution  in  each  position
 *      (columns) and time(rows), in the numpy.loadtxt format. Record 1 per
 *      10 time steps to  avoid  massive  large  files  and  to  force  the
 *      smallness of time step for better precision, so that what is showed
 *      in the results was actually generated with a higher precision.
 *
 * ../gp_data/fileId_domain.dat
 *
 *      text file containning position and time domain information
 *      (1) x1
 *      (2) x2
 *      (3) M
 *      (4) dt
 *      (5) N
 *
 * ========================================================================= */





void SaveConf(char timeinfo, char fnameIn [], char fnameOut [], int Nlines)
{

    int
        i,
        N,
        k,
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
        Vname [30],
        fname [120];

    FILE
        * confFileIn,
        * eqFileIn,
        * confFileOut;

    strcpy(fname, "setup/GP_");
    strcat(fname, fnameIn);
    strcat(fname, "_domain.dat");

    confFileIn = fopen(fname, "r");
    if (confFileIn == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        return;
    }

    strcpy(fname, "setup/GP_");
    strcat(fname, fnameIn);
    strcat(fname, "_eq.dat");

    eqFileIn = fopen(fname, "r");
    if (eqFileIn == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        return;
    }

    strcpy(fname, "../gp_data/");
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

        k = fscanf(confFileIn, "%lf %lf %d %lf %d",
                   &xi, &xf, &Mdx, &dt, &N);

        fprintf(confFileOut, "%d %.15lf %.15lf %.15lf %d ",
                Mdx, xi, xf, dt, N);

        fprintf(confFileOut, "%.15lf %.15lf %.15lf ", a2, imag, g);
        
        fprintf(confFileOut, "%s %.15lf %.15lf %.15lf", Vname, p1, p2, p3);

        fprintf(confFileOut, "\n");
    }

    fclose(eqFileIn); fclose(confFileIn); fclose(confFileOut);

}





void SetupParams(FILE * paramFile, FILE * confFile, Rarray x, double * a2,
     double complex * a1, double * g, Rarray V, double * dt, int * N)
{

/** Read line by line of _domain file and _eq to setup a new integration **/

    int
        k,
        Mdx;

    char
        Vname[30];

    double
        p1,
        p2,
        p3,
        xi,
        xf,
        dx,
        imag;



    // Setup spatial and time domain
    // -----------------------------

    k = fscanf(confFile, "%lf %lf %d %lf %d", &xi, &xf, &Mdx, dt, N);

    // Setup Equation parameters
    // -------------------------

    k = fscanf(paramFile, "%lf %lf %lf %s %lf %lf %lf",
        a2, &imag, g, Vname, &p1, &p2, &p3);

    * a1 = 0 + imag * I;
    
    GetPotential(Mdx + 1, Vname, x, V, p1, p2, p3);
    
    printf("\n\nConfiguration Done\n");

}





int main(int argc, char * argv[])
{

    /*  DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER ARCHITECTURE
     *  --------------------------------------------------------------- */

    mkl_set_num_threads(omp_get_max_threads() / 2);
    omp_set_num_threads(omp_get_max_threads() / 2);



    /* ==================================================================== *
     *
     *                         VARIABLES DEFINITION
     *
     * ==================================================================== */



    int
        i,
        j,
        M,      // # of intervals in spacial domain (sizeof(x) - 1)
        N,      // # of time steps to evolve
        Nlines, // # of initial data to evolve
        trash,  // returned values from scanf (unused)
        cyclic, // boolean-like to set boundary conditions
        method; // method of integrator



    double
        start,      // start trigger to measure time
        time_used,  // end of section with time being measured
        x1,         // lower bound of spacial domain
        x2,         // upper bound of spacial domain
        dx,         // spacial step
        dt,         // time step
        a2,         // Coefficient of (d^2 / d x^2)
        inter,      // interaction strength
        real,       // read real data from file
        imag,       // read imaginary data from file
        * V,        // Potential computed at discretized positions
        * x;        // Array of discretized positions



    double complex
        a1;         // Coefficient of (d / dx)



    char
        timeinfo,
        fname[100], // name to look for files with parameters
        strnum[30]; // string of # of execution



    FILE
        * domain_file,
        * orb_file,
        * eq_file,
        * E_file;



    Carray
        S, // Starts with initial solution, ends with final time-step
        E; // Energy of the system on each time-step





    // Check if there is the right number of command line arguments
    // ------------------------------------------------------------
    if (argc < 5 || argc > 7)
    {
        printf("\nInvalid number of command line arguments, ");
        printf("expected at least 5 and at most 6.\n\n");
        return -1;
    }
    // ------------------------------------------------------------





    // Read data from command line arguments
    // ----------------------------------------------------------------

    timeinfo = argv[1][0]; // character i for imaginary or r for real time

    if (timeinfo != 'r' && timeinfo != 'R')
    {
        if (timeinfo != 'i' && timeinfo != 'I')
        {
            printf("\n\n\tInvalid first argument.\n");
            return -1;
        }
    }

    sscanf(argv[2], "%d", &cyclic);
    if (argc > 5) { sscanf(argv[5], "%d", &method); }
    else          { method = 1;                     }
    if (argc > 6) { sscanf(argv[6], "%d", &Nlines); }
    else          { Nlines = 1;                     }
    // ----------------------------------------------------------------





    /*  ===============================================================
     
               SETUP VALUES FROM FILES - FIXED FOR ALL EXECUTIONS
     
        ===============================================================  */





    // search and read file with values of domain
    // ----------------------------------------------------------------
    strcpy(fname, "setup/GP_");
    strcat(fname, argv[3]);
    strcat(fname, "_domain.dat");

    printf("\nLooking for %s", fname);

    domain_file = fopen(fname, "r");

    if (domain_file == NULL)
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    } else
    {
        printf(" ... Found !\n");
    }

    trash = fscanf(domain_file, "%lf %lf %d", &x1, &x2, &M);

    fclose(domain_file);
    // ----------------------------------------------------------------
    
    
    
    // Setup discretized positions and potential
    // ----------------------------------------------------------------
    dx = (x2 - x1) / M;
    x  = rarrDef(M + 1);
    rarrFillInc(M + 1, x1, dx, x);
    V  = rarrDef(M + 1);
    S  = carrDef(M + 1); // solution
    // ----------------------------------------------------------------



    // write configuration parameters to files
    // ----------------------------------------------------------------
    SaveConf(timeinfo, argv[3], argv[4], Nlines);










    /*  ===============================================================
     
                    LET FILES OPEN TO EXECUTION LIST OF JOBS
     
        ===============================================================  */



    // open file with values of equation parameters
    // ----------------------------------------------------------------
    strcpy(fname, "setup/GP_");
    strcat(fname, argv[3]);
    strcat(fname, "_eq.dat");

    printf("\nLooking for %s", fname);

    eq_file = fopen(fname, "r");

    if (eq_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }
    else
    {
        printf(" ....... Found !\n");
    }
    // ----------------------------------------------------------------



    // open file with values of domain
    // ----------------------------------------------------------------
    strcpy(fname, "setup/GP_");
    strcat(fname, argv[3]);
    strcat(fname, "_domain.dat");

    domain_file = fopen(fname, "r");

    if (domain_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }
    // ----------------------------------------------------------------



    // open file with values of initial condition
    // ----------------------------------------------------------------
    strcpy(fname, "setup/GP_");
    strcat(fname, argv[3]);
    strcat(fname, "_init.dat");

    printf("\nLooking for %s", fname);

    orb_file = fopen(fname, "r");

    if (orb_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }
    else
    {
        printf(" ..... Found !\n");
    }
    
    // open file to write energy values
    // ----------------------------------------------------------------
    strcpy(fname, "../gp_data/");
    strcat(fname, argv[4]);
    strcat(fname, "_energy_imagtime.dat");

    E_file = fopen(fname, "w");

    if (E_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }



    SetupParams(eq_file, domain_file, x, &a2, &a1, &inter, V, &dt, &N);










    /*  ===============================================================
     
                       SETUP INITIAL DATA TO PROPAGATE
     
        ===============================================================  */

    E = carrDef(N + 1); // energy to record progress of convergence

    for (i = 0; i < M + 1; i++)
    {
        trash = fscanf(orb_file, " (%lf%lfj)", &real, &imag);
        S[i] = real + I * imag;
    }
    
    if (timeinfo == 'i' || timeinfo == 'I') fclose(orb_file);

    printf("\nGot Initial condition. Calling time evolution routine\n");










    /*  ===============================================================
     
                             CALL INTEGRATION ROUTINE
     
        ===============================================================  */



    start = omp_get_wtime();


    if (timeinfo == 'r' || timeinfo == 'R')
    {
        SepLine();
        printf("\nDoing real time integration  #%d\n\n", 1);

        strcpy(fname, "../gp_data/");
        strcat(fname, argv[4]);
        strcat(fname, "_line-1_function_realtime.dat");
        switch (method)
        {
            case 1:
                GPCNSMRK4(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                fname, 20);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear CN-SM linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 2:
                GPFFTRK4(M + 1, N, dx, dt, a2, a1, inter, V, S, fname, 20);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear / FFT linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 3:
                GPCNSM(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                fname, 20);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-SM)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 4:
                GPCNLU(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                fname, 20);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-LU)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 5:
                GPFFT(M + 1, N, dx, dt, a2, a1, inter, V, S, fname, 20);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(FFT)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 6:
                CFDS(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                fname, 20);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(CFDS)");
                printf(" : %.3f seconds\n", time_used);
                break;
        }
    }


    if (timeinfo == 'i' || timeinfo == 'I')
    {
        SepLine();
        printf("\nDoing imaginary time integration #%d\n\n", 1);
        switch (method)
        {
            case 1:
                N = IGPCNSMRK4(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear/CN-SM linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 2:
                N = IGPFFTRK4(M + 1, N, dx, dt, a2, a1, inter, V, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear/FFT linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 3:
                N = IGPCNSM(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-SM)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 4:
                N = IGPCNLU(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-LU)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 5:
                N = IGPFFT(M + 1, N, dx, dt, a2, a1, inter, V, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(FFT)");
                printf(" : %.3f seconds\n", time_used);
                break;
        }
    }





    /*  ===============================================================
     
                                   RECORD DATA
     
        ===============================================================  */



    printf("\nRecording data ...\n");

    if (timeinfo == 'i' || timeinfo == 'I')
    {
        strcpy(fname, "../gp_data/");
        strcat(fname, argv[4]);
        strcat(fname, "_line-1_orb_imagtime.dat");

        carr_txt(fname, M + 1, S);
        
        /*
        strcpy(fname, "../gp_data/");
        strcat(fname, argv[4]);
        strcat(fname, "_line-1_E_imagtime.dat");

        carr_txt(fname, N, E);
        */
        fprintf(E_file, "%.10E\n", creal(E[N-1]));
    }





    for ( i = 1; i < Nlines; i ++)
    {

        printf("\n\n\n\n\n");

        // number of line reading in _conf.dat and _eq.dat files
        sprintf(strnum, "%d", i + 1);

        // free old solution
        free(E);

        // read new parameters(one more line) to do another time propagation
        SetupParams(eq_file, domain_file, x, &a2, &a1, &inter, V, &dt, &N);

        E = carrDef(N + 1); // energy to record progress of convergence

        if (timeinfo == 'r' || timeinfo == 'R')
        {

        //  replace initial condition reading more  M + 1 discretized
        //  positions function values from file, for real time domain

            for (j = 0; j < M + 1; j++)
            {
                trash = fscanf(orb_file, " (%lf%lfj)", &real, &imag);
                S[j] = real + I * imag;
            }
        }





        start = omp_get_wtime();

        if (timeinfo == 'r' || timeinfo == 'R')
        {
            SepLine();
            printf("\nDoing real time integration  #%d \n\n", i + 1);

            strcpy(fname, "../gp_data/");
            strcat(fname, argv[4]);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_orb_realtime.dat");

            switch (method)
            {
                case 1:
                    GPCNSMRK4(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                    fname, 10);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear CN-SM linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 2:
                    GPFFTRK4(M + 1, N, dx, dt, a2, a1, inter, V, S,
                    fname, 10);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear / FFT linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 3:
                    GPCNSM(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                    fname, 10);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-SM)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 4:
                    GPCNLU(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                    fname, 10);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-LU)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 5:
                    GPFFT(M + 1, N, dx, dt, a2, a1, inter, V, S, fname, 10);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(FFT)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 6:
                    CFDS(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S,
                    fname, 20);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(CFDS)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
            }
        }


        if (timeinfo == 'i' || timeinfo == 'I')
        {
            SepLine();
            printf("\nDoing imaginary time integration  #%d\n\n", i + 1);
            switch (method)
            {
                case 1:
                    N = IGPCNSMRK4(M + 1, N, dx, dt, a2, a1, inter, V,
                        cyclic, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear/CN-SM linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 2:
                    N = IGPFFTRK4(M + 1, N, dx, dt, a2, a1, inter, V, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear/FFT linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 3:
                    N = IGPCNSM(M + 1, N, dx, dt, a2, a1, inter, V,
                        cyclic, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-SM)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 4:
                    N = IGPCNLU(M + 1, N, dx, dt, a2, a1, inter, V,
                        cyclic, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-LU)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 5:
                    N = IGPFFT(M + 1, N, dx, dt, a2, a1, inter, V, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(FFT)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
            }
        }
        
        printf("\nRecording data ...\n");

        if (timeinfo == 'i' || timeinfo == 'I')
        {
            strcpy(fname, "../gp_data/");
            strcat(fname, argv[4]);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_orb_imagtime.dat");

            carr_txt(fname, M + 1, S);

            /*
            strcpy(fname, "../gp_data/");
            strcat(fname, argv[4]);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_E_imagtime.dat");

            carr_txt(fname, N, E);
            */
            
            fprintf(E_file, "%.10E\n", creal(E[N-1]));
        }
    }





    /* release memory
     * ------------------------------------------------------------------- */
    fclose(eq_file);
    fclose(E_file);
    if (timeinfo == 'r' || timeinfo == 'R') fclose(orb_file);
    fclose(domain_file);
    free(x);
    free(V);
    free(S);
    free(E);
    /* ------------------------------------------------------------------- */



    /*   ==========================    END    ==========================   */

    printf("\n\n");
    return 0;
}
