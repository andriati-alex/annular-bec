#ifndef _NewtonCG_h
#define _NewtonCG_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "calculus.h"
#include "array_memory.h"
#include "iterative_solver.h"
#include "rmatrix_operations.h"

void ncg(int M, double dx, double a2, double complex a1, double inter,
         double mu, Rarray V, Carray f0, Carray S, int * ni, int * cgi);

#endif
