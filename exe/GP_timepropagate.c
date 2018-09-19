#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../include/GP_realtime_integrator.h"
#include "../include/GP_imagtime_integrator.h"





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





int main(int argc, char * argv[])
{

    /*  DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER
     *  -------------------------------------------------- */

    mkl_set_num_threads(omp_get_max_threads() / 2);
    omp_set_num_threads(omp_get_max_threads() / 2);



    /* ==================================================================== *
    /*
     *                         VARIABLES DEFINITION
     *
     * ==================================================================== */

    int
        i,
        M,      // # of intervals in spacial domain (sizeof(x) - 1)
        N,      // # of time steps to evolve
        trash,  // returned values from scanf (unused)
        cyclic, // boolean-like to set boundary conditions
        timeId, // what to do: real ( == 0) or imag ( != 0) time
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
        lambda,     // Strength of delta barrier
        real,       // read real data from file
        imag,       // read imaginary data from file
        * V,        // Potential computed at discretized positions
        * x;        // Array of discretized positions

    double complex
        a1;         // Coefficient of (d / dx)

    char
        fname_in[60],   // name to look for files with parameters
        fname_out[60];  // name of files to write data obtained

    FILE
        * eq_setup_file,
        * out_data;

    Cmatrix S; // To store solution at each time and position step





    /* Check if there is the right number of command line arguments
     * ------------------------------------------------------------ */
    if (argc < 4 || argc > 6)
    {
        printf("\nInvalid number of command line arguments, ");
        printf("expected at least 3 and at most 4.\n\n");
        return -1;
    }
    /* ------------------------------------------------------------ */



    /*  ===============================================================  *
     *
     *                     SETUP VALUES TO VARIABLES
     *
     *  ===============================================================  */



    /* search and read file with values of domain
     * ----------------------------------------------------------------  */
    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_domain.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    trash = fscanf(eq_setup_file, "%lf %lf %d", &x1, &x2, &M);

    fclose(eq_setup_file);
    /* ----------------------------------------------------------------  */



    /* Read data from command line arguments
     * ----------------------------------------------------------------  */
    sscanf(argv[1], "%lf", &dt);
    sscanf(argv[2], "%d",  &N);
    if (argc == 6)
    {
        sscanf(argv[4], "%d", &method);
        sscanf(argv[5], "%d", &timeId);
    }
    else
    {
        timeId = 1;
        if (argc == 5) { sscanf(argv[4], "%d", &method); }
        else           { method = 1;                     }
    }
    /* ----------------------------------------------------------------  */



    /* Setup discretized positions and potential
     * ----------------------------------------------------------------  */
    dx = (x2 - x1) / M;
    x  = rarrDef(M + 1);
    rarrFillInc(M + 1, x1, dx, x);
    V  = rarrDef(M + 1);
    /* ----------------------------------------------------------------  */



    /* search and read file with values of equation parameters
     * ----------------------------------------------------------------  */
    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_eq.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    trash = fscanf(eq_setup_file, "%lf %lf %lf %lf %d",
                   &a2, &imag, &inter, &lambda, &cyclic);

    fclose(eq_setup_file);

    a1 = 0 + imag * I;

    rarrFill(M + 1, 0, V);
    V[M/2] = lambda / dx;  // implement with delta barrier case lambda != 0

    printf("\nEquation coef. and domain successfully setted up.\n");
    /* ----------------------------------------------------------------  */



    /* Read from file and setup initial condition
     * ----------------------------------------------------------------  */
    S = cmatDef(N + 1, M + 1); // matrix to store time step solution

    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_init.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    for (i = 0; i < M + 1; i++)
    {
        trash = fscanf(eq_setup_file, " (%lf+%lfj)", &real, &imag);
        S[0][i] = real + I * imag;
    }

    fclose(eq_setup_file);

    printf("\nGot Initial condition. Calling time evolution routine ...\n");
    /* ----------------------------------------------------------------  */





    /*  ===============================================================  *
     *
     *                        CALL INTEGRATOR ROUTINE
     *
     *  ===============================================================  */

    start = omp_get_wtime();

    if (timeId > 0)
    {
        printf("\n\n\t=================================================\n\n");
        printf("\t * Doing real time integration.\n\n");
        switch (method)
        {
            case 1:
                GPCNSM_all(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-SM)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 2:
                GPCNLU_all(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-LU)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 3:
                GPCNSMRK4_all(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-RK4)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 4:
                GPFFT_all(M + 1, N, dx, dt, a2, a1, inter, V, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(FFT)");
                printf(" : %.3f seconds\n", time_used);
                break;
        }
    }
    else
    {
        printf("\n\n\t=================================================\n\n");
        printf("\t * Doing imaginary time integration.\n\n");
        switch (method)
        {
            case 1:
                IGPCNSM_all(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-SM)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 2:
                IGPCNLU_all(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-LU)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 3:
                IGPCNSMRK4_all(M + 1, N, dx, dt, a2, a1, inter, V, cyclic, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(Crank-Nicolson-RK4)");
                printf(" : %.3f seconds\n", time_used);
                break;
            case 4:
                GPFFT_all(M + 1, N, dx, dt, a2, a1, inter, V, S);
                time_used = (double) (omp_get_wtime() - start);
                printf("\nTime taken to solve(FFT)");
                printf(" : %.3f seconds\n", time_used);
                break;
        }
    }





    /*  ===============================================================  *
     *
     *                             RECORD DATA
     *
     *  ===============================================================  */

    strcpy(fname_out, "../gp_data/");
    strcat(fname_out, argv[3]);
    strcat(fname_out, "_time.dat");

    printf("\nRecording data ...\n");
    cmat_txt(fname_out, N + 1, 10, M + 1, 1, S);

    strcpy(fname_out, "../gp_data/");
    strcat(fname_out, argv[3]);
    strcat(fname_out, "_domain.dat");

    out_data = fopen(fname_out, "w");

    if (out_data == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_out); return -1; }

    fprintf(out_data, "%.15lf %.15lf %d %.10lf %d", x1, x2, M, dt, N);

    fclose(out_data);




    /* release memory
     * ------------------------------------------------------------------- */
    free(x);
    free(V);
    cmatFree(N + 1, S);
    /* ------------------------------------------------------------------- */



    /*   ==========================    END    ==========================   */

    printf("\n\n");
    return 0;
}
