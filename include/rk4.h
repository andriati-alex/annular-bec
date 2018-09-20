#ifndef _rk4_h
#define _rk4_h

#include "array_memory.h"

/* DESCRIPTION
 *
 * Evolve from an initial condition x one time step in x_step with fourth
 * order Runge-Kutta routine                                           */

void RK4step(int M, double dt, double t, Carray x, Carray extra, Carray x_step,
     void (*dxdt)(int, double , Carray, Carray, Carray));

#endif
