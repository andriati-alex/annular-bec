#ifndef _rk4_h
#define _rk4_h

#include "array_memory.h"
#include "matrix_operations.h"

/* DESCRIPTION
 *
 * Evolve from an initial condition y one time step in y_step with fourth
 * order Runge-Kutta routine, with linear dependecy on function  over the
 * vector y, the function calls  then  must be described by matrix vector
 * multiplication, whereas the matrix is passed in CCS format in F.    */
void RK4step(int M, double dt, CCSmat F, Carray y, Carray y_step);


/* Call the routine RK4step N times. Setup the solution in Matrix S */
void RK4(int M, int N, double dt, CCSmat F, Cmatrix S);

#endif
