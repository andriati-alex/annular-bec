#ifndef _MCTDHB_integrator_h
#define _MCTDHB_integrator_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "array_memory.h"
#include "matrix_operations.h"
#include "tridiagonal_solver.h"
#include "coef_routines.h"
#include "calculus.h"



/* ====================================================================== */
/*                                                                        */
/*                          DATA-TYPE DEFINITION                          */
/*                                                                        */
/* ====================================================================== */



struct _MCTDHBsetup
{
    int 
        Mpos,   // # of discretized positions (# divisions + 1)
        Morb,   // # of orbitals
        Npar,   // # of particles
        ** IF;  // IF[i] point to the occupation number vetor of C[i]
    long
        nc,         // Total # of configurations of Fock states
        ** NCmat;   // NCmat[n][m] # with n particles / m orbitals
    double
        dx,     // space step
        xi,     // initial position discretized value
        xf,     // final position discretized value
        a2,     // factor multiplying d2 / dx2
        inter,  // know as g, contact interaction strength
        * V;    // Array with the values of one-particle potential
    double complex
        a1;     // factor multiplying d / dx (pure imaginary)
};

typedef struct _MCTDHBsetup * MCTDHBsetup;



/* ====================================================================== */
/*                                                                        */
/*                          FUNCTION PROTOTYPES                           */
/*                                                                        */
/* ====================================================================== */



MCTDHBsetup AllocMCTDHBdata(
        int Npar, 
        int Morb,
        int Mpos,
        double xi,
        double xf,
        double a2,
        double inter,
        double * V,
        double complex a1);

void EraseMCTDHBdata(MCTDHBsetup MC);

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

void RK4step(MCTDHBsetup MC, Cmatrix orb, Carray C, double dt);

void RHSforRK4(MCTDHBsetup MC, Carray C, Cmatrix Orb,
               Cmatrix Ho, Carray Hint,
               Carray newC, Cmatrix newOrb);

void LinearPart(MCTDHBsetup MC, CCSmat rhs_mat, Carray upper, Carray lower,
                Carray mid, Cmatrix Orb);

void MCTDHB_time_evolution(MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt,
        int cyclic);

#endif
