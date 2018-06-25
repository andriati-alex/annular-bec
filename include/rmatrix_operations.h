#ifndef _rmatrix_operations_h
#define _rmatrix_operations_h

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "array.h"
#include "array_memory.h"



/*          ***********************************************          */
/*                   SETUP VALUES IN MATRIX ENTRIES                  */
/*          ***********************************************          */



void rmatFill(int m, int n, double z, Rmatrix M);
/* Fill all matrix with a constant value
 * *************************************
 *
 *  m the number of lines
 *  n the number of columns
 *  z the value to fill the entire matrix
 *
 * *************************************/

void rmatFillDK(int n, int k, Rarray z, Rmatrix M);
/* Fill k Diagonal of square matrix of n by n elements with an array z
 * *******************************************************************
 *
 * z must contain n - k elements
 * k = 0 is the main diagonal
 * k > 0 upper diagonals
 * k < 0 lower diagonals
 *
 * *******************************************************************/

void rmatFillTri(int n, Rarray upper, Rarray mid, Rarray lower, Rmatrix M);
/* Fill a matrix with just the tridiagonals entries, rest with zeros */

RCCSmat triToRCCS(int n, Rarray upper, Rarray lower, Rarray mid);
/* Fill a Sparse Matrix in Compressed Column Storage format from a tridiagonal
 * ***************************************************************************
 *
 * n the dimension of the matrix. Returns a pointer to CCS struct
 *
 *****************************************************************************/

RCCSmat CyclicToRCCS(int n, Rarray upper, Rarray lower, Rarray mid);
/* Fill in CCS format given tridiagonal cyclic system */

void setValueRCCS(int n, int i, int j, int col, double z, RCCSmat M);
/* Set a value in the i row and original column number col */

void addValueRCCS(int n, int i, int j, double z, RCCSmat M);
/* Add a value in the i row and CCS column(not the true index) position j */

RCCSmat emptyRCCS(int n, int max_nonzeros);
/* Initialize an empty CCS matrix */


/*          **********************************************          */
/*          MATRIX-VECTOR AND MATRIX-MATRIX MULTIPLICATION          */
/*          **********************************************          */



void rmatvec(int m, int n, Rmatrix M, Rarray v, Rarray ans);
/* General Matrix Vector multiplication: M . v = ans
 * *************************************************
 *
 *  m number of lines of M
 *  n number of columns of M and components of v
 *
 * *************************************************/

void RCCSvec(int n, Rarray vals, int * restrict cols, int m, Rarray vec,
            Rarray ans);
/* Matrix(in CCS format) vector multiplication
 * 
 * Given CCSmat A the arguments taken are
 *
 * vals = A->vec
 * cols = A->col
 * m    = A->m
 *
 * *******************************************/

void rmatmat(int m, int n, int l, Rmatrix M, Rmatrix A, Rmatrix ans);
/* General Matrix Matrix multiplication: M . A = ans
 * *************************************************
 *
 *  M has m(lines) by n(columns)
 *  A has n(lines) by l(columns)
 *  ans has m(lines) by l(columns)
 *
 * *************************************************/

#endif
