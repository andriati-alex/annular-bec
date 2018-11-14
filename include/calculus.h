#ifndef _calculus_h
#define _calculus_h

#include <mkl.h>
#include <mkl_dfti.h>
#include <math.h>

#include "array_memory.h"
#include "array_operations.h"

void Ortonormalize(int Mfun, int Mpos, double dx, Cmatrix F);

// Renormalize f such that || f || = norm
void renormalize(int n, Carray f, double dx, double norm);

/* DERIVATIVES AND INTEGRATION ROUTINES */

double complex Csimps(int n, Carray f, double dx);

double Rsimps(int n, Rarray f, double dx);

void dxFFT(int n, Carray f, double dx, Carray dfdx);
/* Automatically periodic boundary due to use of FFT */

void dxCyclic(int n, Carray f, double dx, Carray dfdx);
/* Compute using finite centered differences with periodic boundary */

double complex Functional(int M, double dx, double a2, double complex a1,
                          double inter, Rarray V, Carray f);

void applyL0(int n, Carray f, double dx, double a2, double complex a1, 
             Rarray V, double inter, double mu, Carray L0f);

double complex GPkinect(int M, double a2, double complex a1, double dx,
               Carray psi);

double complex GPtrap(int M, Rarray V, double dx, Carray psi);

double GPinter(int M, double g, double dx, Carray psi);

double complex GPvirial(int M, double a2, double complex a1, double g,
               Rarray V, double dx, Carray psi);

#endif
