#ifndef _calculus_h
#define _calculus_h

#include <mkl.h>
#include <mkl_dfti.h>
#include <math.h>

#include "array_memory.h"
#include "array_operations.h"

#define PI 3.141592653589793 // to compute derivatives with FFT

/* DERIVATIVES AND INTEGRATION ROUTINES */

double complex Csimps(int n, Carray f, double dx);

double Rsimps(int n, Rarray f, double dx);

void dxFFT(int n, Carray f, double dx, Carray dfdx);
/* Automatically periodic boundary due to use of FFT */

void dxCyclic(int n, Carray f, double dx, Carray dfdx);
/* Compute using finite centered differences with periodic boundary */

#endif