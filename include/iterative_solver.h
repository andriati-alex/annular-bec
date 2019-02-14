#ifndef _iterative_solver_h
#define _iterative_solver_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "array_memory.h"
#include "array_operations.h"
#include "matrix_operations.h"
#include "tridiagonal_solver.h"



int CCG(int n, CCSmat A, Carray b, Carray x, double eps, int maxiter, 
               Carray upper, Carray lower, Carray mid);
/* As in MpreCCG solve tridiagonal system to apply pre-conditioning */



int RCG(int n, RCCSmat A, Rarray b, Rarray x, double eps, int maxiter,
        Rarray upper, Rarray lower, Rarray mid);
/* Use pre-conditioned CG but for real entries */



#endif
