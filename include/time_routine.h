#ifndef _time_routine_h
#define _time_routine_h

#include <mkl.h>
#include <mkl_dfti.h>

#include "tridiagonal_solver.h"
#include "matrix_operations.h"
#include "array_operations.h"

#define PI 3.141592653589793 // for Fourier frequencies



/*                         ***********************                         */
/*                         TIME PROPAGATE ROUTINES                         */
/*                         ***********************                         */

/* GENERAL DECRIPTION
 *
 * Use split-step scheme to separate non-linearity  and treat the linear
 * part in diffetent methods, Crank-Nicolson  and  spectral (using FFT),
 * whereas the first apply Sherman-Morrison or LU decomposition to solve
 * the linear system required.
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


void CNsm(int M, int N, double dx, double dt, double a2, double complex a1,
          double inter, Rarray V, int cyclic, Cmatrix S);
/* Crank-Nicolson with Sherman-Morrison for linear part */

void CNlu(int M, int N, double dx, double dt, double a2, double complex a1,
          double inter, Rarray V, int cyclic, Cmatrix S);
/* Crank-Nicolson with modified LU decomposition */

void spectral(int M, int N, double dx, double dt, double a2, double complex a1,
              double inter, Rarray V, Cmatrix S);
/* Use MKL for fourier tranform and compute derivatives */

void lastStep(int M, int N, double dx, double dt, double a2, double complex a1,
              double inter, Rarray V, int cyclic, Carray S);
/* Return just an array with the last time step obtained */

void ApplyRHS(int M, double inter, Carray S, Carray rhs);

void RK4gp(int M, int N, double dx, double dt, double a2, double complex a1,
           double inter, Rarray V, int cyclic, Cmatrix S);

#endif
