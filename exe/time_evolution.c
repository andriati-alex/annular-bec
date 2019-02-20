#include <string.h>
#include <stdio.h>
#include <math.h>
#include "realtime_integrator.h"
#include "imagtime_integrator.h"
#include "linear_potential.h"





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





void ReachNewLine(FILE * f)
{

    // Read until get new line in a opened file.

    char
        sentinel;

    while (1)
    {
        fscanf(f, "%c", &sentinel);
        if (sentinel == '\n' || sentinel == EOF) return;
    }
}





void SaveConf(FILE * f, EqDataPkg EQ, double dt, double N)
{

    double imag;

    imag = cimag(EQ->a1);

    fprintf(f, "%d %.15lf %.15lf %.10lf %d ", EQ->Mpos, EQ->xi, EQ->xf, dt, N);

    fprintf(f, "%.15lf %.15lf %.15lf ", EQ->a2, imag, EQ->inter);

    fprintf(f, "%.15lf %.15lf %.15lf", EQ->p[0], EQ->p[1], EQ->p[2]);

    fprintf(f, "\n");

}





EqDataPkg SetupParams(FILE * paramFile, FILE * confFile,
          char Vname[], double * dt, int * N)
{

/** Read line by line of _domain file and _eq to setup a new integration **/

    int
        k,
        Mdx;

    double
        g,
        a2,
        xi,
        xf,
        imag,
        p[3];



    // Read spatial and time domain settings
    // -------------------------------------

    k = fscanf(confFile, "%lf %lf %d %lf %d", &xi, &xf, &Mdx, dt, N);

    // Read a line of numbers corresponding to equation parameters
    // -----------------------------------------------------------

    k = fscanf(paramFile,"%lf %lf %lf %lf %lf %lf", &a2, &imag, &g,
        &p[0], &p[1], &p[2]);

    return PackEqData(Mdx + 1, xi, xf, a2, g, I * imag, Vname, p);

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
        rate,
        M,      // # of intervals in spacial domain (sizeof(x) - 1)
        N,      // # of time steps to evolve
        Nlines, // # of initial data to evolve
        trash,  // returned values from scanf (unused)
        cyclic, // boolean-like to set boundary conditions
        method, // method of integrator
        resetinit;



    double
        start,      // start trigger to measure time
        time_used,  // end of section with time being measured
        dt,         // time step
        real,       // read real data from file
        imag;       // read imaginary data from file



    char
        c,
        timeinfo,
        fname[100], // name to look for files with parameters
        strnum[30], // string of # of execution
        potname[50],
        infname[100],
        outfname[100];



    FILE
        * domain_file,
        * job_file,
        * orb_file,
        * eq_file,
        * E_file;



    Carray
        S, // Starts with initial solution, ends with final time-step
        E; // Energy of the system on each time-step



    EqDataPkg EQ;





    job_file = fopen("job.conf", "r");

    if (job_file == NULL) // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n", "job.conf");
        exit(EXIT_FAILURE);
    }

    i = 1;

    while ( (c  = getc(job_file)) != EOF)
    {

        // jump comment line
        if (c == '#') { ReachNewLine(job_file); continue; }
        else          { fseek(job_file, -1, SEEK_CUR);    }

        switch (i)
        {
            case 1:
                fscanf(job_file, "%s", fname);
                timeinfo = fname[0];
                i = i + 1;
                break;
            case 2:
                fscanf(job_file, "%s", potname);
                i = i + 1;
                break;
            case 3:
                fscanf(job_file, "%d", &cyclic);
                i = i + 1;
                break;
            case 4:
                fscanf(job_file, "%s", infname);
                i = i + 1;
                break;
            case 5:
                fscanf(job_file, "%s", outfname);
                i = i + 1;
                break;
            case 6:
                fscanf(job_file, "%d", &method);
                i = i + 1;
                break;
            case 7:
                fscanf(job_file, "%d", &Nlines);
                i = i + 1;
                break;
            case 8:
                fscanf(job_file, "%d", &resetinit);
                i = i + 1;
                break;
        }

        ReachNewLine(job_file);

    }

    fclose(job_file);





    if (timeinfo != 'r' && timeinfo != 'R')
    {
        if (timeinfo != 'i' && timeinfo != 'I')
        {
            printf("\n\n\tInvalid type of propagation.\n");
            exit(EXIT_FAILURE);
        }
    }










    /*  ===============================================================
     
                   LET FILES OPENNED TO EXECUTE LIST OF JOBS
     
        ===============================================================  */



    // open file with values of equation parameters
    // ----------------------------------------------------------------
    strcpy(fname, "input/");
    strcat(fname, infname);
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
    strcpy(fname, "input/");
    strcat(fname, infname);
    strcat(fname, "_domain.dat");

    printf("\nLooking for %s", fname);

    domain_file = fopen(fname, "r");

    if (domain_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }
    else
    {
        printf(" ... Found !\n");
    }
    // ----------------------------------------------------------------



    // open file with values of initial condition
    // ----------------------------------------------------------------
    strcpy(fname, "input/");
    strcat(fname, infname);
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
    strcpy(fname, "output/");
    strcat(fname, outfname);
    strcat(fname, "_energy.dat");

    E_file = fopen(fname, "w");

    if (E_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }



    // open file to write parameters of domain and equation
    // ----------------------------------------------------------------
    strcpy(fname, "output/");
    strcat(fname, outfname);

    if (timeinfo == 'i' || timeinfo == 'I')
    {
        strcat(fname, "_conf_imagtime.dat");
    } else
    {
        strcat(fname, "_conf_realtime.dat");
    }

    job_file = fopen(fname, "w");

    if (job_file == NULL)  // impossible to open file
    {
        printf("\n\n\tERROR: impossible to open file %s\n\n", fname);
        return -1;
    }

    fprintf(job_file, "# Trap Id : %s\n", potname);



    EQ = SetupParams(eq_file, domain_file, potname, &dt, &N);
    M  = EQ->Mpos - 1;
    S  = carrDef(M + 1); // solution at discretized positions

    rate = N / 1000 + 1;










    /*  ===============================================================
     
                       SETUP INITIAL DATA TO PROPAGATE
     
        ===============================================================  */

    E = carrDef(N + 1); // energy at each time step

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
        sepline();
        printf("\nDoing real time integration  #%d\n\n", 1);

        strcpy(fname, "output/");
        strcat(fname, outfname);
        strcat(fname, "_line-1_orb_realtime.dat");
        switch (method)
        {
            case 1:
                SSCNRK4(EQ, N, dt, cyclic, S, fname, rate);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear CN-SM linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 2:
                SSFFTRK4(EQ, N, dt, S, fname, rate);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear / FFT linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 3:
                SSCNSM(EQ, N, dt, cyclic, S, fname, rate);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-SM)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 4:
                SSCNLU(EQ, N, dt, cyclic, S, fname, rate);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-LU)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 5:
                SSFFT(EQ, N, dt, S, fname, rate);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(FFT)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 6:
                CFDS(EQ, N, dt, cyclic, S, fname, rate);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(CFDS)");
                printf(" : %.3f seconds\n", time_used);
                break;
        }
    }


    if (timeinfo == 'i' || timeinfo == 'I')
    {
        sepline();
        printf("\nDoing imaginary time integration #%d\n\n", 1);
        switch (method)
        {
            case 1:
                N = ISSCNRK4(EQ, N, dt, cyclic, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear/CN-SM linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 2:
                N = ISSFFTRK4(EQ, N, dt, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(RK4 nonlinear/FFT linear)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 3:
                N = ISSCNSM(EQ, N, dt, cyclic, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-SM)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 4:
                N = ISSCNLU(EQ, N, dt, cyclic, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-LU)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 5:
                N = ISSFFT(EQ, N, dt, S, E);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(FFT)");
                printf(" : %.3f seconds\n", time_used);
                break;
        }

        // Record data
        strcpy(fname, "output/");
        strcat(fname, outfname);
        strcat(fname, "_line-1_orb_imagtime.dat");

        carr_txt(fname, M + 1, S);

        fprintf(E_file, "%.10E\n", creal(E[N-1]));
    }

    SaveConf(job_file, EQ, dt, N);





    for ( i = 1; i < Nlines; i ++)
    {

        printf("\n\n\n\n\n\n\n\n");
        
        ReleaseEqDataPkg(EQ);

        free(E);

        // number of line reading in _conf.dat and _eq.dat files
        sprintf(strnum, "%d", i + 1);

        // read new parameters(one more line) to do another time propagation
        EQ = SetupParams(eq_file, domain_file, potname, &dt, &N);

        rate = N / 1000 + 1;

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
        else
        {

            // Set the same initial condition for all executions if
            // resetinit is True for imaginary time propagation

            if (resetinit)
            {

                strcpy(fname, "input/");
                strcat(fname, infname);
                strcat(fname, "_init.dat");

                printf("\nUsing the same initial condition");

                orb_file = fopen(fname, "r");

                for (j = 0; j < M + 1; j++)
                {
                    trash = fscanf(orb_file, " (%lf%lfj)", &real, &imag);
                    S[j] = real + I * imag;
                }

                fclose(orb_file);
            }
        }





        start = omp_get_wtime();

        if (timeinfo == 'r' || timeinfo == 'R')
        {
            sepline();
            printf("\nDoing real time integration  #%d \n\n", i + 1);

            strcpy(fname, "output/");
            strcat(fname, outfname);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_orb_realtime.dat");

            switch (method)
            {
                case 1:
                    SSCNRK4(EQ, N, dt, cyclic, S, fname, rate);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear CN-SM linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 2:
                    SSFFTRK4(EQ, N, dt, S, fname, rate);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear / FFT linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 3:
                    SSCNSM(EQ, N, dt, cyclic, S, fname, rate);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-SM)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 4:
                    SSCNLU(EQ, N, dt, cyclic, S, fname, rate);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-LU)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 5:
                    SSFFT(EQ, N, dt, S, fname, rate);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(FFT)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 6:
                    CFDS(EQ, N, dt, cyclic, S, fname, rate);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(CFDS)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
            }
        }


        if (timeinfo == 'i' || timeinfo == 'I')
        {
            sepline();
            printf("\nDoing imaginary time integration  #%d\n\n", i + 1);
            switch (method)
            {
                case 1:
                    N = ISSCNRK4(EQ, N, dt, cyclic, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear/CN linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 2:
                    N = ISSFFTRK4(EQ, N, dt, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(RK4 nonlinear/FFT linear)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 3:
                    N = ISSCNSM(EQ, N, dt, cyclic, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-SM)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 4:
                    N = ISSCNLU(EQ, N, dt, cyclic, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(Crank-Nicolson-LU)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
                case 5:
                    N = ISSFFT(EQ, N, dt, S, E);
                    time_used = (double) (omp_get_wtime() - start);
                    printf("\nTime taken to solve(FFT)");
                    printf(" : %.3f seconds\n", time_used);
                    break;
            }

            strcpy(fname, "output/");
            strcat(fname, outfname);
            strcat(fname, "_line-");
            strcat(fname, strnum);
            strcat(fname, "_orb_imagtime.dat");

            carr_txt(fname, M + 1, S);

            fprintf(E_file, "%.10E\n", creal(E[N-1]));

        }

        SaveConf(job_file, EQ, dt, N);
    }





    /* release memory
     * ------------------------------------------------------------------- */
    fclose(job_file);
    fclose(eq_file);
    fclose(E_file);
    if (timeinfo == 'r' || timeinfo == 'R') fclose(orb_file);
    fclose(domain_file);
    free(S);
    free(E);
    ReleaseEqDataPkg(EQ);
    /* ------------------------------------------------------------------- */



    /*   ==========================    END    ==========================   */

    printf("\n\n");
    return 0;
}
