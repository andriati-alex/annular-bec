#ifndef _matrix_operations_h
#define _matrix_operations_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "array.h"
#include "array_memory.h"
#include "array_operations.h"



/*          ***********************************************          */
/*                   SETUP VALUES IN MATRIX ENTRIES                  */
/*          ***********************************************          */



void cmatFill(int m, int n, double complex z, Cmatrix M);
/* Fill all matrix with a constant value
 * *************************************
 *
 *  m the number of lines
 *  n the number of columns
 *  z the value to fill the entire matrix
 *
 * *************************************/



void cmatFillDK(int n, int k, Carray z, Cmatrix M);
/* Fill k Diagonal of square matrix of n by n elements with an array z
 * *******************************************************************
 *
 * z must contain n - k elements
 * k = 0 is the main diagonal
 * k > 0 upper diagonals
 * k < 0 lower diagonals
 *
 * *******************************************************************/



void cmatFillTri(int n, Carray upper, Carray mid, Carray lower, Cmatrix M);
/*   Fill a matrix with just the tridiagonals entries, rest with zeros   */



CCSmat triToCCS(int n, Carray upper, Carray lower, Carray mid);
/* Fill a Sparse Matrix in Compressed Column Storage format from a tridiagonal
 * ***************************************************************************
 *
 * n the dimension of the matrix. Returns a pointer to CCS struct
 *
 *****************************************************************************/



CCSmat CyclicToCCS(int n, Carray upper, Carray lower, Carray mid);
/*      Fill in CCS format given tridiagonal cyclic system      */



void setValueCCS(int n, int i, int j, int col, double complex z, CCSmat M);
/*        Set a value in the i row and original column number col        */



CCSmat emptyCCS(int n, int max_nonzeros);
/*    Initialize an empty CCS matrix   */



void RowMajor(int m, int n, Cmatrix M, Carray v);
/* Store Matrix M(m x n) in a vector v using Row Major scheme */

CCSmat CNmat(int M, double dx, double complex dt, double a2,
     double complex a1, double inter, Rarray V, int cyclic,
     Carray upper, Carray lower, Carray mid);



/*          **********************************************          */
/*          MATRIX-VECTOR AND MATRIX-MATRIX MULTIPLICATION          */
/*          **********************************************          */



void cmatvec(int m, int n, Cmatrix M, Carray v, Carray ans);
/* General Matrix Vector multiplication: M . v = ans
 * *************************************************
 *
 *  m number of lines of M
 *  n number of columns of M and components of v
 *
 * *************************************************/



void CCSvec(int n, Carray vals, int * restrict cols, int m, Carray vec,
            Carray ans);
/* Matrix(in CCS format) vector multiplication
 * 
 * Given CCSmat A the arguments taken are
 *
 * vals = A->vec
 * cols = A->col
 * m    = A->m
 *
 * *******************************************/



void cmatmat(int m, int n, int l, Cmatrix M, Cmatrix A, Cmatrix ans);
/* General Matrix Matrix multiplication: M . A = ans
 * *************************************************
 *
 *  M has m(lines) by n(columns)
 *  A has n(lines) by l(columns)
 *  ans has m(lines) by l(columns)
 *
 * *************************************************/

int HermitianInv(int M, Cmatrix A, Cmatrix A_inv);
/* Invert an hermitian matrix */



/*          ***********************************************          */
/*          DETERMINANT AND INVERSION OF TRIDIAGONAL SYSTEM          */
/*          ***********************************************          */



double complex detTrik(int n, Carray upper, Carray lower, Carray mid);
/* Use recursion to compute determinant of block submatrix of tridiagonal
 *
 * Starting from the first row and column take k rows and columns to
 * define another square matrix and compute its determinant
 *
 * **********************************************************************/



Carray detTri(int n, Carray upper, Carray lower, Carray mid);
/* Take a block submatrix k by k and store its determinant in k-component */



Carray phiTri(int n, Carray upper, Carray lower, Carray mid);
/* Backwards determinant. Useful to compute inverse of tridiagonal matrix */



void invTri(int n, Carray upper, Carray lower, Carray mid, Cmatrix ans);
/* Compute the inverse of tridiagonal matrix
 * *****************************************
 *
 * detTri stores all theta components from 
 * determinant recurrence phiTri is a kind 
 * of backward determinant recurrence.
 * Search tridiagonal matrices in Wikipedia 
 * for implementation details.
 *
 * *****************************************/

#endif
