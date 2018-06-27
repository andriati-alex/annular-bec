#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../include/NewtonCG.h"

/* 
 * OBTAIN AN STEADY STATE USING NEWTON METHOD FOR OPERATORS
 * ********************************************************
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
 *      (2) imaginary part of first order derivative
 *      (3) interaction strength
 *      (4) potential strength
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
 * mu fileId
 *
 *      mu     -> Value of chemical potential
 *      fileId -> the prefix of required file name
 *
 * CALL
 * ****
 *
 * ./time_evolution mu fileId
 *
 * OUTPUT FILES
 * ************
 *
 * **************************************************************************/

int main(int argc, char * argv[])
{

    /* DEFINE THE NUMBER OF THREADS BASED ON THE COMPUTER */
    mkl_set_num_threads(omp_get_max_threads() / 2);
    omp_set_num_threads(omp_get_max_threads() / 2);

    double start, time_used; // show time taken for time routine solver

    int trash; // Useless returned values

    if (argc != 3)
    {
        printf("\nInvalid number of command line arguments, ");
        printf("expected 2.\n\n");
        return -1;
    }



    /*          ********************************************          */
    /*          Setup domain of solution and method to solve          */
    /*          ********************************************          */



    double x1, x2, dx;  // Position domain = [x1, x2] interval
    double mu;          // Chemical potential
    unsigned int M;     // M is number of dx, N number of time steps
    Rarray x;           // Discretized positions Vector

    char fname_in[40];      // file configuration names
    FILE * eq_setup_file;   // pointer to file

    /*** Command line arguments ***/

    sscanf(argv[1], "%lf", &mu);

    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[2]);
    strcat(fname_in, "_domain.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    trash = fscanf(eq_setup_file, "%lf %lf %d", &x1, &x2, &M);

    fclose(eq_setup_file); // finish reading of file

    dx = (x2 - x1) / M;
    x  = rarrDef(M + 1);
    rarrFillInc(M + 1, x1, dx, x);



    /*                    *************************                    */
    /*                    Setup equation parameters                    */
    /*                    *************************                    */
    


    double a2,                  // second order derivative coef.
           inter,               // interaction strength coef.
           lambda,              // potential parameter
           a1imag;

    double complex a1;          // First order derivative coef.

    Rarray V = rarrDef(M + 1);  // Potential in discretized positions
    
    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[2]);
    strcat(fname_in, "_eq.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    trash = fscanf(eq_setup_file, "%lf %lf %lf %lf",
                   &a2, &a1imag, &inter, &lambda);

    fclose(eq_setup_file); // finish reading of file

    a1 = 0 + a1imag * I;

    rarrFill(M + 1, 0, V);
    V[M/2] = lambda / dx;

    printf("\nEquation coef. and domain successfully setted up.\n");



    /*              *************************************              */
    /*              Evolve solution from an initial state              */
    /*              *************************************              */



    int ni, cgi; // Count number of iterations
    double real, imag; // to read data from file

    Carray f0 = carrDef(M + 1); // Initial attempt
    Carray S  = carrDef(M + 1); // Store solution

    strcpy(fname_in, "setup/");
    strcat(fname_in, argv[2]);
    strcat(fname_in, "_init.dat");

    printf("\nLooking for %s\n", fname_in);

    eq_setup_file = fopen(fname_in, "r");

    if (eq_setup_file == NULL)  // impossible to open file
    { printf("ERROR: impossible to open file %s\n", fname_in); return -1; }

    for (int i = 0; i < M + 1; i++)
    {
        trash = fscanf(eq_setup_file, " (%lf+%lfj)", &real, &imag);
        f0[i] = real + I * imag;
    }

    fclose(eq_setup_file); // finish the reading of file

    printf("\nGot Initial attempt. Calling NewtonCG routines ...\n");

    start  = omp_get_wtime();
    time_used = (double) (omp_get_wtime() - start);
    ncg(M, dx, a2, a1, inter, mu, V, f0, S, &ni, &cgi);
    printf("\nTime taken NewtonCG : %.3f seconds\n", time_used);



    /*                          ***********                          */
    /*                          Record data                          */
    /*                          ***********                          */



    /*** release memory ***/

    free(x); free(V); free(f0); free(S);

    /* END */

    printf("\n");
    return 0;
}
