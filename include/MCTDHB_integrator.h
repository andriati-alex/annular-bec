#ifndef _MCTDHB_integrator_h
#define _MCTDHB_integrator_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "array_memory.h"
#include "matrix_operations.h"
#include "coef_routines.h"
#include "calculus.h"



/* ====================================================================== */
/*                                                                        */
/*                          DATA-TYPE DEFINITION                          */
/*                                                                        */
/* ====================================================================== */



struct orbitals
{
    int Morb;   // # of orbitals
    int Mpos;   // # discretized positions (# of divisions + 1)
    double dx;  // spatial step
    double xi;
    double xf;
    Cmatrix O;  // Values of orbitals on each position (Morb x M)
};

typedef struct orbitals * MCorbital;

struct coeficients
{
    int Npar;       // # of particles
    int Morb;       // # of orbitals
    long nc;        // all possible fock states
    int ** IF;      // All occu. vector by row
    long ** NCmat;  // All outcomes from combinatorial
    Carray C;       // Array of coeficients
};

typedef struct coeficients * MCcoef;

struct eq_parameters
{
    double a2;          // Second order derivative term
    double inter;       // interaction 'strength'
    double complex a1;  // First order derivative term
    Rarray V;           // Values of potential in each position
};

typedef struct eq_parameters * EqSetup;



/* ====================================================================== */
/*                                                                        */
/*                          FUNCTION PROTOTYPES                           */
/*                                                                        */
/* ====================================================================== */



MCcoef AllocCoef(int Npar, int Morb);
/* **********************************
 *
 * Return the pointer to the struct of Coeficients
 * -----------------------------------------------
 * 
 * Configure all fields in the struct but C.
 *
 * ***********************************************/

EqSetup AllocEq(double a2, double complex a1, double inter, Rarray V);
/* *******************************************************************
 *
 * Return the pointer to the struct of Equation setup parameters
 * -------------------------------------------------------------
 * 
 * a1 multiplies (d / dx)
 * a2 multiplies (d2 / dx2)
 * inter multiplies the nonlinearity
 * V have the values of single-particle potential at each position
 * 
 * *******************************************************************/

void SetupHo(int Morb, int Mpos, Cmatrix Omat, double dx,
             double a2, double complex a1, Rarray V, Cmatrix Ho);
/* **************************************************************
 *
 * Setup matrix elements of one-body part on Ho
 * --------------------------------------------
 * 
 * **************************************************************/

void SetupHint(int Morb, int Mpos, Cmatrix Omat, double dx,
               double inter, Carray Hint);
/* ********************************************************
 *
 * Setup matrix elements of two-body part on Hint
 * ----------------------------------------------
 * 
 * ********************************************************/

double complex Proj_Hint(int M, int k, int i,
                         Cmatrix rho_inv, Carray rho2, Carray Hint);

double complex NonLinear(int M, int k, int n,
                         Cmatrix Omat, Cmatrix rho_inv, Carray rho2);

void RK4step(MCorbital psi, MCcoef coef, EqSetup eq, double dt);

void RHSforRK4(int N, int M, long ** NCmat, int ** IF,
               int Mpos, double dx, EqSetup eq,
               Carray C, Cmatrix Orb,
               Cmatrix Ho, Carray Hint,
               Carray newC, Cmatrix newOrb);

#endif
