#ifndef _iterative_solver_h
#define _iterative_solver_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "array_memory.h"
#include "array_operations.h"
#include "matrix_operations.h"
#include "rmatrix_operations.h"
#include "tridiagonal_solver.h"



int CCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps);
/* Solve (complex)linear system by Conjugate-Gradient method
 * *********************************************************
 *
 * n is the size of the system
 * A symmetric matrix (self-adjoint)
 * A . x = b (RHS)
 * x initial guess ends up with solution
 * eps the tolerated residual error
 * ---------------------------------------------------------
 *
 * RETURN the number of iterations to converge
 *
 * *********************************************************/



int preCCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps,
           Cmatrix M);
/* Pre-Conditioned Conjugate-Gradient Method
 * *****************************************
 *
 * takes an extra(last) argument that is 
 * properly the inversion of some Matrix M
 * similar to A, but with known inversion.
 * Here the method needs an extra matrix 
 * vector multiplication but can iterate 
 * less times in return.
 *
 * *****************************************/



int MpreCCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps, 
            Carray upper, Carray lower, Carray mid);
/* Instead of take the inversion of M, takes M itself as tridiagonal mat
 * *********************************************************************
 *
 * As the inversion is not given, in the loop it solves tridiagonal system
 *
 * *********************************************************************/



int CCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps);
/* Take A as Compressed Column Storaged(CCS) matrix */



int preCCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps,
              Cmatrix M);
/* As in preCCG use matrix vector multiplication by pre-conditioning */



int MpreCCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps, 
               Carray upper, Carray lower, Carray mid);
/* As in MpreCCG solve tridiagonal system to apply pre-conditioning */



int RCG(unsigned int n, RCCSmat A, Rarray b, Rarray x, double eps,
        Rarray upper, Rarray lower, Rarray mid);
/* Use pre-conditioned CG but for real entries */

#endif
