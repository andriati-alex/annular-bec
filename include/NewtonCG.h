#ifndef _NewtonCG_h
#define _NewtonCG_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <limits.h>
#include "GP_functional.h"
#include "iterative_solver.h"
#include "rmatrix_operations.h"

struct IterNCG{
    int newton;
	double cg_mean;
    double cg_std;
    int cg_total;
    int cg_max;
    int cg_min;
};

void iteration_info(struct IterNCG It);

struct IterNCG ncg(int M, double tol, int maxiter, double dx, double a2,
                   double complex a1, double inter, double mu, Rarray V,
                   Carray f0);

#endif
