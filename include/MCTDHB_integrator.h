#ifndef _MCTDHB_integrator_h
#define _MCTDHB_integrator_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <string.h>

#include "MCTDHB_datatype.h"
#include "MCTDHB_observables.h"
#include "MCTDHB_configurations.h"
#include "matrix_operations.h"
#include "tridiagonal_solver.h"



/* ======================================================================== */
/*                                                                          */
/*                           FUNCTION PROTOTYPES                            */
/*                                                                          */
/* ======================================================================== */





void applyHconf(MCTDHBsetup MC, Carray C, Cmatrix Ho, Carray Hint, Carray out);
/* Give the state coefficients of a state (out) after apply the many-body
 * Hamiltonian on a state whose  coeficients  in  the  occupation  number
 * basis are C[i]. Ho contains the matrix elements of  one-body  operator
 * and Hint contains the matrix elements of two-body operator.        **/





/* ==========================================================================
 *                                                                   *
 *          Function to apply nonlinear part of obitals PDE          *
 *          -----------------------------------------------          *
 *                                                                   */





double complex nonlinear (int M, int k, int n, double g, Cmatrix Orb,
               Cmatrix Rinv, Carray R2, Cmatrix Ho, Carray Hint );





/* ==========================================================================
 *                                                                   *
 *                    Time Evolution of equations                    *
 *                    ---------------------------                    *
 *                                                                   */





void OrbDDT (MCTDHBsetup MC, Carray C, Cmatrix Orb, Cmatrix newOrb,
     Cmatrix ho, Carray Hint);





void OrbConfDDT
(   // The Right-Hand-Side(RHS) of a system of Differential equations
    MCTDHBsetup MC,
    Carray C,
    Cmatrix Orb,
    Cmatrix Ho,
    Carray Hint,
    Carray newC,   // Values after the operations have been applied
    Cmatrix newOrb // Values after the operations have been applied
);





void lanczos(MCTDHBsetup MCdata, Cmatrix Ho, Carray Hint,
     int lm, Carray diag, Carray offdiag, Cmatrix lvec);





void RK4orbstep (MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt);
void ABMorb (MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt);
void lanczosCstep (MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt);





void RK4lanczosBefore
(   // Evolve nonlinear part of orbitals coupled with coeficients
    MCTDHBsetup MC,
    Cmatrix orb, // End up modified by the evolution
    Carray C,    // End up modified by the evolution
    double dt
);

void RK4lanczosAfter
(   // Evolve nonlinear part of orbitals coupled with coeficients
    MCTDHBsetup MC,
    Cmatrix orb, // End up modified by the evolution
    Carray C,    // End up modified by the evolution
    double dt
);





void IRK4step
(   // Evolve nonlinear part with imaginary time
    MCTDHBsetup MC,
    Cmatrix Orb,
    Carray C,
    double complex dt
);





void RK4step
(   // Evolve nonlinear part
    MCTDHBsetup MC,
    Cmatrix Orb,
    Carray C,
    double dt
);





void LinearPartSM
(   // Evolve a time-step the linear part of PDE(orbitals)
    int Mpos,
    int Morb,
    CCSmat cnmat, // Matrix from CN approach(discretization)
    Carray upper, // Upper diagonal
    Carray lower, // Lower diagonal
    Carray mid,   // Main diagonal of tridiagonal system
    Cmatrix Orb
);





void LinearPartLU (int Mpos, int Morb, CCSmat cnmat, Carray upper,
     Carray lower, Carray mid, Cmatrix Orb);





/* ==========================================================================
 *                                                                   *
 *                     Main routine to be called                     *
 *                    ---------------------------                    *
 *                                                                   */





void MCTDHB_CN_REAL_ABM (MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt,
     int Nsteps, int method, int cyclic, char fname [], int n);

void MCTDHB_CN_REAL (MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt,
     int Nsteps, int method, int cyclic, char fname[], int n);





void MCTDHB_CN_IMAG (MCTDHBsetup MC, Cmatrix Orb, Carray C, Carray E,
     double dT, int Nsteps, int cyclic);






#endif
