#ifndef _GP_realtime_integrator_h
#define _GP_realtime_integrator_h

#include <mkl.h>
#include <mkl_dfti.h>

#include "tridiagonal_solver.h"
#include "matrix_operations.h"
#include "array_operations.h"





/*  ======================================================================  *
 *
 *           REAL TIME INTEGRATOR ROUTINES FOR GROSS-PITAEVSKII
 *
 *  ======================================================================  */





/* =======================================================================
 *
 * GENERAL DESCRIPTION
 * -------------------
 *
 * Use split-step scheme to separate non-linearity  and treat the linear
 * part in diffetent methods, Crank-Nicolson  and  spectral  (using FFT),
 * whereas the first apply Sherman-Morrison or LU decomposition to solve
 * the linear system required.
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
 *
 * ======================================================================= */

void record_step(FILE * f, int M, Carray v);



void GPCNSM(int M, int N, double dx, double dt, double a2,
     double complex a1, double inter, Rarray V, int cyclic, Carray S);
/* Crank-Nicolson with Sherman-Morrison for linear part */





void GPCNLU(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname[], int n);
/* Crank-Nicolson with modified LU decomposition */





void GPFFT(int M, int N, double dx, double dt, double a2,
     double complex a1, double inter, Rarray V, Carray S);
/* Use MKL for fourier tranform and compute derivatives */





void DDT(int M, double inter, Carray S, Carray rhs);
/* the right-hand-side of derivatives from nonlinear part */





void GPCNSMRK4_all(int M, int N, double dx, double dt, double a2,
     double complex a1, double inter, Rarray V, int cyclic, Cmatrix S);

#endif
