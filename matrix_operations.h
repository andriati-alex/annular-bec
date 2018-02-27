#ifndef _matrix_operations_h
#define _matrix_operations_h

#include "array.h"
#include "array_memory.h"

void cmatFill(unsigned int m, unsigned int n, double complex z, Cmatrix M);
/* Fill all matrix with a constant value
 * *************************************
 *
 *  m the number of lines
 *  n the number of columns
 *  z the value to fill the entire matrix
 *
 ***************************************/

void cmatFillDK(unsigned int n, int k, Carray z, Cmatrix M);
/* Fill k Diagonal of square matrix of n by n elements with an array z
 *********************************************************************
 *
 * z must contain n - k elements
 * k = 0 is the main diagonal
 * k > 0 upper diagonals
 * k < 0 lower diagonals
 *
 *********************************************************************/

void cmatFillTri(unsigned int n, 
                 Carray upper, Carray mid, Carray lower, Cmatrix M);
/* Tridiagonal Matrix
 ********************
 *
 * n size of the main diagonal
 * upper has size n - 1
 * lower has size n - 1
 * Other positions of M are filled with zero
 *
 ********************/

struct CCS *triToCCS(unsigned int n, Carray upper, Carray lower, Carray mid);
/* Fille a Sparse Matrix in Compressed Column Storage format from a tridiagonal
 ******************************************************************************
 *
 * n the dimension of the matrix
 * returns a struct containing the vector of values, column index and
 * a integer that is the maximum number of nonzero elements in a row.
 *
 *****************************************************************************/

struct CCS *CyclicToCCS(unsigned int n, Carray upper, Carray lower, Carray mid);
/* Fill a Sparse Matrix in Compressed Column Storage format from a Cyclic */

void CCSvec(unsigned int n, Carray vals, int * cols, int m, 
            Carray vec, Carray ans);
/* Matrix vector multiplication(in CCS format) 
 * *******************************************
 *
 * ans must be filled with zeros
 *
 * *******************************************/

void cmatvec(unsigned int m, unsigned int n, Cmatrix M, Carray v, Carray ans);
/* Matrix Vector multiplication: M . v = ans
 * *****************************************
 *
 *  m number of lines of M
 *  n number of columns of M and components of v
 *
 *******************************************/

void cmatmat(unsigned int m, unsigned int n, unsigned int l,
             Cmatrix M, Cmatrix A, Cmatrix ans);
/* Matrix Matrix multiplication: M . A = ans
 *******************************************
 *
 *  M has m(lines) by n(columns) dimension
 *  A has n(lines) by l(columns)
 *  ans has m(lines) by l(columns)
 *
 *******************************************/

void cmatvecTri(unsigned int n, Cmatrix M, Carray v, Carray ans);
/* tridiagonal-Matrix vector multiplication */

void cmatvecDiag(unsigned int n, int *U, int *L, Cmatrix M,
                 Carray v, Carray ans);
/* Multiply a Matrix with only diagonals filled by a vector
 **********************************************************
 *
 * U is vector with first line column's position of upper diagonals
 * U[i] >=0 ( = 0 is interpreted as the main diagonal.
 * U[i] < 0 stop looking for diagonals.
 * Example U = { 0, 2, 5, -1} has main, second upper and fith upper diags
 *
 * L is vector with first column line's position of lower diagonals
 * L[i] >= 1
 *
 **********************************************************/

double complex detTrik(Carray upper, Carray lower, Carray mid, unsigned int k);
/* With recursion compute determinant of block submatrix of k by k elements */

Carray detTri(Carray upper, Carray lower, Carray mid, unsigned int n);
/* Take a block submatrix k by k and store its determinant in k-component */

Carray phiTri(Carray upper, Carray lower, Carray mid, unsigned int n);
/* Backwards determinant. Useful to compute inverse of tridiagonal matrix */

void invTri(Carray upper, Carray lower, Carray mid,
            unsigned int n, Cmatrix ans);
/* Compute the inverse of tridiagonal matrix
 *******************************************
 *
 * detTri stores all theta components from determinant recurrence
 * phiTri is a kind of backward determinant recurrence
 * Search tridiagonal matrices in Wikipedia for implementation details.
 *
 *******************************************/

#endif
