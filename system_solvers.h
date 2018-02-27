#ifndef _system_solvers_h
#define _ststem_solvers_h

#include "array.h"
#include "array_memory.h"
#include "array_operations.h"
#include "matrix_operations.h"





double cond(unsigned int n, Carray upper, Carray lower, Carray mid);
double errBack(unsigned int n, Carray upper, Carray lower, Carray mid);
/* Compute condition number u
 * **************************
 *
 * See BUENO, M. Isabel & DOPICO, FROILÁN M. for the meaning of condition
 * numbers in "Stability and sensitivity of tridiagonal LU factorization"
 *
 * A good condition number is < 1E10 (assure at least 6 digits precision).
 *
 * The Backward error should be of 10E-15
 *
 ****************************/





double complex checkLU(unsigned int n, Carray upper, Carray lower, Carray mid);
/* Check if the tridiagonal system has L.U. factorization
 * ******************************************************
 *
 * n is the size of the system
 * upper is upper diagonal and has size (n-1)
 * lower is lower diagonal and has size (n-1)
 * mid is the main diagonal of size (n)
 *
 * The process ends up with the determinant, otherwise it gives zero
 * (a source of breakdown to LU factorization)
 *
 ********************************************************/





void triDiag(unsigned int n, Carray upper, Carray lower, Carray mid,
             Carray RHS, Carray ans);
/* Solve Tri-diagonal linear system using A = L . U factorization
 * **************************************************************
 *
 * n is the size of the system
 * upper is upper diagonal and has size (n-1)
 * lower is lower diagonal and has size (n-1)
 * mid is the main diagonal of size (n)
 * RHS (Right Hand side) vector of A . x = RHS
 * ans ends up with solution
 *
 ****************************************************************/





void triCyclicLU(unsigned int n, Carray upper, Carray lower, Carray mid,
                 Carray RHS, Carray ans);
/* Solve Cyclic Tri-diagonal linear system using Modified LU decomposition
 * ***********************************************************************
 *
 * upper is upper diagonal and has size n (last element is the top coupling)
 * lower is lower diagonal and has size n (''   ''   is the bottom coupling)
 * mid is the main diagonal of size n
 * RHS (Right Hand side) vector of A . x = RHS
 * n is the size of the system
 * ans ends up with solution
 * -----------------------------------------------------------------------
 *
 * Implementarion follows (Thiab R.)Taha's paper
 * "Solution of Periodic Tridiagonal Linear Systems on a Hypercube"
 * Method fails if first main diagonal element is zero.
 *
 **************************************************************************/





void triCyclicSM(unsigned int n, Carray upper, Carray lower, Carray mid,
                 Carray RHS, Carray ans);
/* Solve Cyclic Tri-diagonal linear system using Sherman-Morrison
 * **************************************************************
 *
 * upper is upper diagonal and has size n (last element is the top coupling)
 * lower is lower diagonal and has size n (''   ''   is the bottom coupling)
 * mid is the main diagonal of size n
 * RHS (Right Hand side) vector of A . x = RHS
 * n is the size of the system
 * ans ends up with solution
 * --------------------------------------------------------------
 *
 * Implementation follows numerical recipes, pag ~ 70
 *
 *****************************************************************/





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
 ***********************************************************/





int preCCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps,
           Cmatrix M);
/* Pre-Conditioned Conjugate-Gradient Method
 * *****************************************
 *
 * takes an extra(last) argument that is properly the inversion
 * of some Matrix M. Here the method needs an extra matrix vector
 * multiplication
 *
 *******************************************/





int MpreCCG(unsigned int n, Cmatrix A, Carray b, Carray x, double eps, 
            Carray upper, Carray lower, Carray mid);
/* Instead of take the inversion of M, takes M itself as tridiagonal mat
 * *********************************************************************
 *
 * As the inversion is not given, in the loop it solves tridiagonal system
 *
 ***********************************************************************/





int CCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps);
/* Take A as Compressed Column Storaged(CCS) matrix */





int preCCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps,
              Cmatrix M);
/* As in preCCG use matrix vector multiplication by pre-conditioning
 * *****************************************************************
 *
 * Further More uses A as a CCS matrix
 *
 *******************************************************************/





int MpreCCSCCG(unsigned int n, CCSmat A, Carray b, Carray x, double eps, 
            Carray upper, Carray lower, Carray mid);
/* As in MpreCCG solve tridiagonal system to apply pre-conditioning
 * ****************************************************************
 *
 * Further More uses A as a CCS matrix
 *
 ******************************************************************/

#endif
