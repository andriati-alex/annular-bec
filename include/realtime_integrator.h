#ifndef _realtime_integrator_h
#define _realtime_integrator_h

#include <mkl.h>
#include <mkl_dfti.h>

#include "tridiagonal_solver.h"
#include "matrix_operations.h"
#include "array_operations.h"
#include "observables.h"
#include "inout.h"
#include "rk4.h"
#include "data_structure.h"





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





void SSCNSM(EqDataPkg, int N, double dt, int cyclic, Carray S,
     char fname[], int n);
/* ---------------------------------------------------------
 * Crank-Nicolson with Sherman-Morrison to solve linear part
 * --------------------------------------------------------- */





void SSCNLU(EqDataPkg, int N, double dt, int cyclic, Carray S,
     char fname[], int n);
/* ---------------------------------------------------------
 * Crank-Nicolson with LU decomposition to solve linear part
 * --------------------------------------------------------- */





void SSFFT(EqDataPkg, int N, double dt, Carray S, char fname[], int n);
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





void SSCNRK4(EqDataPkg, int N, double dt, int cyclic, Carray S,
     char fname [], int n);
/* ---------------------------------------
 * Crank-Nicolson with Sherman-Morrison to
 * solve linear part and RK4 to  nonlinear
 * --------------------------------------- */





void SSFFTRK4(EqDataPkg, int N, double dt, Carray S, char fname [], int n);
/* -----------------------------------------------------------
 * Use FFT to solve derivative part and RK4 for potential part
 * ----------------------------------------------------------- */





void CFDS(EqDataPkg, int N, double dt, int cyclic, Carray S,
     char fname [], int n);

#endif
