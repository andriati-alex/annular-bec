#ifndef _itime_routine_h
#define _itime_routine_h

#include <mkl.h>
#include <mkl_dfti.h>

#include "tridiagonal_solver.h"
#include "matrix_operations.h"
#include "calculus.h"

#define PI 3.141592653589793 // to compute derivatives from FFT

/*                    *********************************                    */
/*                    IMAGINARY TIME PROPAGATE ROUTINES                    */
/*                    *********************************                    */

/* GENERAL DECRIPTION
 *
 * Use split-step scheme  to separate non-linearity  and treat the linear
 * part in  diffetent methods, Crank-Nicolson  and  spectral (using FFT),
 * whereas the first apply Sherman-Morrison or LU decomposition  to solve
 * the linear system required. In each time step renormalize the solution
 * fixing from initial guess the norm.
 *
 * ARGUMENTS
 *
 * The first 4 arguments are discretization of domain in space and time
 * M size of positions(M - 1 intervals of dx), N number  of  time steps
 * evolved with time-step size dt. Next we have the equation coeficients
 * ak that multiply  the  k-th order derivative, inter multiply the non-
 * linearity,  a  vector with position dependent potential, cyclic as a
 * boolean to put periodic or zero boundary conditions(just for CN), and
 * finally S has position in columns and time steps in rows. Case S is
 * array then return just the N-th time step(the last one).
 *
 */

void iCNSM(int M, int N, double dx, double dT, double a2, double complex a1,
           double inter, Rarray V, int cyclic, Cmatrix S);

void iCNLU(int M, int N, double dx, double dT, double a2, double complex a1,
           double inter, Rarray V, int cyclic, Cmatrix S);

void ispectral(int M, int N, double dx, double dT, double a2,
               double complex a1, double inter, Rarray V, Cmatrix S);

void iSMlaststep(int M, int N, double dx, double dT, double a2,
                 double complex a1, double inter, Rarray V, int cyclic,
                 Carray S);

#endif
