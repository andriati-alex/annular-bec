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



/* ======================================================================== */
/*                                                                          */
/*                           DATA-TYPE DEFINITION                           */
/*                                                                          */
/* ======================================================================== */



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



/* ======================================================================== */
/*                                                                          */
/*                           FUNCTION PROTOTYPES                            */
/*                                                                          */
/* ======================================================================== */



MCTDHBsetup AllocMCTDHBdata
(
    int Npar,
    int Morb,
    int Mpos,
    double xi,
    double xf,
    double a2,
    double inter,
    double * V,
    double complex a1
);

void EraseMCTDHBdata(MCTDHBsetup MC);



/* ==========================================================================
 *                                                                *
 *        Setup matrix elements of one-body Hamiltonian Ho        *
 *        ------------------------------------------------        *
 *                                                                */

void SetupHo
(
    int Morb,
    int Mpos,
    Cmatrix Omat,
    double dx,
    double a2,
    double complex a1,
    Rarray V,
    Cmatrix Ho
);



/* ==========================================================================
 *                                                                  *
 *        Setup matrix elements of two-body Hamiltonian Hint        *
 *        --------------------------------------------------        *
 *                                                                  */

void SetupHint
(
    int Morb,
    int Mpos,
    Cmatrix Omat,
    double dx,
    double inter,
    Carray Hint
);



/* ==========================================================================
 *                                                                   *
 *          Function to apply nonlinear part of obitals PDE          *
 *          -----------------------------------------------          *
 *                                                                   */

double complex Proj_Hint
(   // The one that comes from Projective part
    int M, // Total # of orbitals
    int k, // A fixed orbital k
    int i, // A fixed orbital i
    Cmatrix rho_inv,
    Carray rho2,
    Carray Hint
);

double complex NonLinear
(   // The one that comes from identity part
    int M, // Total # of orbitals
    int k, // A fixed orbital
    int n, // A given discretized position
    Cmatrix Omat,
    Cmatrix rho_inv,
    Carray rho2
);



/* ==========================================================================
 *                                                                   *
 *                    Time Evolution of equations                    *
 *                    ---------------------------                    *
 *                                                                   */

void RK4step
(   // Evolve nonlinear part of orbitals coupled with coeficients
    MCTDHBsetup MC,
    Cmatrix orb, // End up modified by the evolution
    Carray C,    // End up modified by the evolution
    double dt
);

void RHSforRK4
(   // The Right-Hand-Side(RHS) of a system of Differential equations
    MCTDHBsetup MC,
    Carray C,
    Cmatrix Orb,
    Cmatrix Ho,
    Carray Hint,
    Carray newC,   // Values after the operations have been applied
    Cmatrix newOrb // Values after the operations have been applied
);

void LinearPart
(   // Evolve a time-step the linear part of PDE(orbitals)
    MCTDHBsetup MC,
    CCSmat rhs_mat, // Matrix from CN approach(discretization)
    Carray upper,   // Upper diagonal
    Carray lower,   // Lower diagonal
    Carray mid,     // Main diagonal of tridiagonal system
    Cmatrix Orb
);

void MCTDHB_time_evolution
(   // Call the subroutines to solve nonlinear and linear part
    MCTDHBsetup MC,
    Cmatrix Orb, // Modified at each time step with the solution(orbitals)
    Carray C,    // Modified at each time step with the solution(coeficients)
    double dt,
    int Nsteps,
    int cyclic
);

#endif
