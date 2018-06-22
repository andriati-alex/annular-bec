#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../include/time_routine.h"

/* 
 * OBTAIN TIME EVOLUTION SOLUTION OF GROSS-PITAEVSKII EQUATION
 * ***********************************************************
 *
 * REQUIRED FILES
 * **************
 *
 * setup/fileId_init.dat
 *
 *      is a text file with an array of complex numbers of size M + 1
 *      where M is the number of space intervals of the domain.  Must
 *      be formatted according to numpy.savetxt function.
 *
 * setup/fileId_eq.dat
 *
 *      is a text file with values of equation coefficients
 *      in the following order:
 *
 *      (1) second order derivative
 *      (2) real part of first order derivative
 *      (3) imaginary part of first order derivative
 *      (4) interaction strength
 *      (5) potential strength
 *
 * setup/fileId_domain.dat
 *
 *      a text file with position domain information over which  the
 *      fileId_init was generated. The numbers is in following order
 *
 *      (1) x_i
 *      (2) x_f
 *      (3) M the number of steps of size (xf - xi) / M
 *
 * COMMAND LINE ARGUMENTS
 * **********************
 * 
 * dt N fileId solverId(optional)
 *
 *      dt       -> the time step
 *      N        -> number of time steps to evolve
 *      fileId   -> the prefix of required file name
 *      solverId -> Method to solve the linear part(optional)
 *
 * PS: solverId, if given, must be:
 *
 *      1 Crank-Nicolson with Sherman-Morrison
 *      2 Crank-Nicolson with special LU-decomposition
 *      3 Fourier Transforms
 *
 * CALL
 * ****
 *
 * ./time_evolution dt N fileId solverId(optional)
 *
 * OUTPUT FILES
 * ************
 *
 * ../gp_data/fileId_time.dat
 *
 *      a text file with matrix-like format with solution in each  position
 *      (columns) and time(rows), in the numpy.loadtxt like format.  Record
 *      1 every 10 time steps to avoid massive large files and to force the
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
 * **************************************************************************/

int main(int argc, char * argv[])
{
    /* DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER */
    mkl_set_num_threads(omp_get_max_threads() / 2);
    omp_set_num_threads(omp_get_max_threads() / 2);

    // Useless returned values
    int trash;

    if (argc < 4 || argc > 5)
    {
        printf("\nInvalid number of command line arguments, ");
        printf("expected at least 3 and at most 4.\n\n");
        return -1;
    }



    /*          ********************************************          */
    /*          Setup domain of solution and method to solve          */
    /*          ********************************************          */



    double x1, x2, dx, dt;  // Position domain = [x1, x2] interval
    unsigned int M, N;      // M is number of dx, N number of time steps
    unsigned int solverId;  // Specify method to solve linear part
    Rarray x;               // Discretized positions Vector

    char fname_in[40];      // file configuration names
    FILE * eq_setup_file;   // pointer to file

    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_domain.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    trash = fscanf(eq_setup_file, "%lf %lf %d", &x1, &x2, &M);

    fclose(eq_setup_file); // finish reading of file

    sscanf(argv[1], "%lf", &dt);
    sscanf(argv[2], "%d",  &N);
    if (argc == 5) { sscanf(argv[4], "%d", &solverId); }
    else           { solverId = 1; }

    dx = (x2 - x1) / M;
    x = rarrDef(M + 1);
    rarrFillInc(M + 1, x1, dx, x);



    /*                    *************************                    */
    /*                    Setup equation parameters                    */
    /*                    *************************                    */



    double a2,                  // second order derivative coef.
           inter,               // interaction strength coef.
           lambda,              // potential parameter
           a1real,
           a1imag;
    double complex a1;          // First order derivative coef.
    Rarray V = rarrDef(M + 1);  // Potential in discretized positions

    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_eq.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    trash = fscanf(eq_setup_file, "%lf %lf %lf %lf %lf",
            &a2, &a1real, &a1imag, &inter, &lambda);

    fclose(eq_setup_file); // finish reading of file

    a1 = a1real + a1imag * I;

    rarrFill(M + 1, 0, V);
    V[M/2] = lambda / dx;

    printf("\nEquation coef. and domain successfully setted up.\n");



    /*              *************************************              */
    /*              Evolve solution from an initial state              */
    /*              *************************************              */



    double start, time_used; // show time taken for time routine solver
    double real, imag;       // to read data from file

    Cmatrix S = cmatDef(N + 1, M + 1); // matrix to store each step solution

    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[3]);
    strcat(fname_in, "_init.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    for (int i = 0; i < M + 1; i++)
    {
        trash = fscanf(eq_setup_file, " (%lf+%lfj)", &real, &imag);
        S[0][i] = real + I * imag;
    }

    fclose(eq_setup_file); // finish the reading of file

    printf("\nGot Initial condition. Calling time evolution routine ...\n");

    start  = omp_get_wtime();
    switch (solverId) {
        case 1:
            CNsm(M + 1, N, dx, dt, a2, a1, inter, V, 1, S);
            time_used = (double) (omp_get_wtime() - start);
            printf("\nTime taken to solve(SM) : %.3f seconds\n", time_used);
            break;
        case 2:
            CNlu(M + 1, N, dx, dt, a2, a1, inter, V, 1, S);
            time_used = (double) (omp_get_wtime() - start);
            printf("\nTime taken to solve(LU) : %.3f seconds\n", time_used);
            break;
        case 3:
            spectral(M + 1, N, dx, dt, a2, a1, inter, V, S);
            time_used = (double) (omp_get_wtime() - start);
            printf("\nTime taken to solve(FFT) : %.3f seconds\n", time_used);
            break;
    }



    /*                          ***********                          */
    /*                          Record data                          */
    /*                          ***********                          */



    char fname_out[30];

    strcpy(fname_out, "../gp_data/");
    strcat(fname_out, argv[3]);
    strcat(fname_out, "_time.dat");

    printf("\nRecording data ...\n");
    cmat_txt(fname_out, N + 1, 10, M + 1, 1, S);

    strcpy(fname_out, "../gp_data/");
    strcat(fname_out, argv[3]);
    strcat(fname_out, "_domain.dat");

    FILE * out_data = fopen(fname_out, "w");

    if (out_data == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_out); return -1; }

    fprintf(out_data, "%.2lf %.2lf %d %.8lf %d", x1, x2, M, dt, N);

    fclose(out_data);

    /* release memory */
    free(x);
    free(V);
    cmatFree(N + 1, S);

    /* END */

    printf("\n");
    return 0;
}
