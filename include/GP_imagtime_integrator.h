#ifndef _GP_imagtime_integrator_h
#define _GP_imagtime_integrator_h

#include <mkl.h>
#include <mkl_dfti.h>

#include "tridiagonal_solver.h"
#include "matrix_operations.h"
#include "calculus.h"





/*                    *********************************                    */
/*                    IMAGINARY TIME PROPAGATE ROUTINES                    */
/*                    *********************************                    */





/* =========================================================================
 *
 * GENERAL DECRIPTION
 * ------------------
 *
 * Use split-step scheme  to separate non-linearity  and treat the linear
 * part in  diffetent methods, Crank-Nicolson  and  spectral  (using FFT),
 * whereas the first apply Sherman-Morrison or LU decomposition  to solve
 * the linear system required. In each time step renormalize the solution
 * fixing the norm from initial guess.
 *
 *
 *
 *
 *
 * ARGUMENTS
 * ---------
 *
 * The first 4 arguments are discretization of domain in space and  time
 * M size of positions(M - 1 intervals of dx), N number  of  time  steps
 * evolved with time-step size dt. Next we have the equation coeficients
 * ak that multiply  the  k-th order derivative, inter multiply the non-
 * linearity,  a  vector with position dependent potential, cyclic  as a
 * boolean to put periodic or zero boundary conditions(just for CN), and
 * finally S has position in columns and time steps in rows. Case  S  is
 * array then return just the N-th time step(the last one).
 *
 * ========================================================================= */





void IGPCNSM(int M, int N, double dx, double dT, double a2, double complex a1,
      double inter, Rarray V, int cyclic, Carray S, Carray E);





void IGPCNLU(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, Carray E);





void IGPFFT(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, Carray S, Carray E);

#endif
