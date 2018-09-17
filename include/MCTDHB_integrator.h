#ifndef _MCTDHB_integrator_h
#define _MCTDHB_integrator_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "MCTDHB_datatype.h"
#include "MCTDHB_observables.h"
#include "MCTDHB_configurations.h"
#include "array_memory.h"
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





void RK4lanczos
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
    CCSmat rhs_mat, // Matrix from CN approach(discretization)
    Carray upper,   // Upper diagonal
    Carray lower,   // Lower diagonal
    Carray mid,     // Main diagonal of tridiagonal system
    Cmatrix Orb
);





/* ==========================================================================
 *                                                                   *
 *                     Main routine to be called                     *
 *                    ---------------------------                    *
 *                                                                   */





void MCTDHB_REAL_LanczosRK4I
(   // Call the subroutines to solve nonlinear and linear part
    MCTDHBsetup MC,
    Cmatrix Orb, // Modified/Overwritten at each time step
    Carray C,    // Modified/Overwritten at each time step
    Carray E,    // Energy at each time step
    double dt,
    int Nsteps,
    int cyclic
);





void MCTDHB_REAL_RK4I
(   // Call the subroutines to solve nonlinear and linear part
    MCTDHBsetup MC,
    Cmatrix Orb, // Modified/Overwritten at each time step
    Carray C,    // Modified/Overwritten at each time step
    Carray E,    // Energy at each time step
    double dt,
    int Nsteps,
    int cyclic
);





void MCTDHB_IMAG_RK4I
(   // Call subroutines with imaginary time to get ground state
    MCTDHBsetup MC,
    Cmatrix Orb, // Modified/Overwritten and renormalized at each time step
    Carray C,    // Modified/Overwritten and renormalized at each time step
    Carray E,    // Energy at each time step
    double dT,
    int Nsteps,
    int cyclic
);






#endif
