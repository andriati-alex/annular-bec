#ifndef _GP_realtime_integrator_h
#define _GP_realtime_integrator_h

#include <mkl.h>
#include <mkl_dfti.h>

#include "tridiagonal_solver.h"
#include "matrix_operations.h"
#include "array_operations.h"
#include "GP_functional.h"
#include "rk4.h"





/*  ======================================================================  *
 *
 *           REAL TIME INTEGRATOR ROUTINES FOR GROSS-PITAEVSKII
 *
 *  ======================================================================  */





/* =======================================================================
 *
 *
 * MODULE OF FUNCTIONS THAT INTEGRATE THE GROSS-PITAEVSKII EQUATION
 *
 *
 * i dS/dt = ( a2 (d^2/dx^2) + a1 (d/dx) + V(x) + inter |S|^2 ) S(x,t)
 *
 *
 * This module implement the mostly well know and used  integrators  to
 * the Gross-Pitaevskii(GP) equation. The routines are named as follows
 *
 *  - First 2 letters : GP (from Gross-Pitaevskii)
 *
 *  - Next letters : identify the methods
 *      (1) CN for Crank-Nicolson finite differences scheme
 *      (2) FFT for Fast Fourier transform to deal with derivatives
 *
 *  - The next letters apply to CN on how to solve the linear system
 *      (1) SM for Sherman-Morrison formula.
 *      (2) LU for the decomposition.
 *  
 *  - For those who use 4th order Runge-Kutta method to integrate the
 *    nonlinear part from the split-step, we have additionally RK4 as
 *    a suffix, otherwise,  it  is  done going one step forward using
 *    integration by rectangles and doing again the same steps  using
 *    integration by trapezium rule
 *
 *  M is the number of discretized points (size of arrays)
 *  N is the number of time-steps to be propagated the initial condition
 *
 *  CN methods supports both cyclic and zero boundary condition as
 *  identified by the cyclic(boolean) parameter.
 *
 * ======================================================================= */



void record_step(FILE * f, int M, Carray v);
/* ------------------------------------------------------------
 * Record v array data (of size M) in a row of a openned file f
 * ------------------------------------------------------------ */





void GPCNSM(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname[], int n);
/* ---------------------------------------------------------
 * Crank-Nicolson with Sherman-Morrison to solve linear part
 * --------------------------------------------------------- */





void GPCNLU(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname[], int n);
/* ---------------------------------------------------------
 * Crank-Nicolson with LU decomposition to solve linear part
 * --------------------------------------------------------- */





void GPFFT(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, Carray S, char fname[], int n);
/* -------------------------------------------------------
 * Use MKL fourier tranform routine to compute derivatives
 * ------------------------------------------------------- */





void NonLinearDDT(int M, double t, Carray Psi, Carray inter, Carray Dpsi);
/* --------------------------------------------------------------------
 * Time-derivative from nonlinear part after split-step (called in RK4)
 * -------------------------------------------------------------------- */





void NonLinearVDDT(int M, double t, Carray Psi, Carray FullPot, Carray Dpsi);
/* ------------------------------------------------------------------------
 * Time-derivative from nonderivative part after split-step (called in RK4)
 * ------------------------------------------------------------------------ */





void GPCNSMRK4(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname [], int n);
/* ---------------------------------------
 * Crank-Nicolson with Sherman-Morrison to
 * solve linear part and RK4 to  nonlinear
 * --------------------------------------- */





void GPFFTRK4(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, Carray S, char fname[], int n);
/* -----------------------------------------------------------
 * Use FFT to solve derivative part and RK4 for potential part
 * ----------------------------------------------------------- */

#endif
