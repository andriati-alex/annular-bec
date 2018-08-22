#ifndef _array_h
#define _array_h

#include <complex.h>
#include <mkl.h>

/* DATATYPES OF THE PACKAGE
 * *********************************************************
 *
 * Description
 * ===========
 *
 * Header with all datatypes used on this project
 *
 *
 *
 * Important Details
 * =================
 *
 * Pointers declared with restrict keyword suppose no 
 * aliasing between pointers. For this you can't pass
 * the same pointer in two different arguments of the 
 * routines.  Moreover this procedure is suitable for
 * parallelized codes.
 *
 * For example don't call: carrAdd(n, v, u, v) aiming
 * produce the effect of v[i] = u[i] + v[i].
 *
 *
 *
 * The Compressed-Column-Storage(CCS) of a Matrix
 * ==============================================
 *
 * This method consist in save memory when  storing a
 * sparse Matrix by shifting all the nonzero elements
 * to the left. Then the resulting matrix has only  m
 * columns,  m being  the maximum  number of  nonzero
 * elements in a same row, and the same number of row
 * than the former one.  Thereafter  this  matrix  is
 * organized in a  vector  following  a  Column Major
 * ordering as well as its former column positions are
 * stored in a vector of integers of the same size.
 *
 * ********************************************************/



/************************ Constants ***********************/

#define LAPACK_ROW_MAJOR 101
#define PI 3.141592653589793



/********************** MKL Pointers **********************/

typedef MKL_Complex16 * CMKLarray;



/***************** Real vector and matrix *****************/

typedef double *  restrict Rarray;
typedef double ** restrict Rmatrix;



/**************** Complex vector and matrix ***************/

typedef double complex *  restrict Carray;
typedef double complex ** restrict Cmatrix;



/******* Compressed Column storage of an n x n Matrix *****/

struct CCS{
    int  m; // max number of non-zero elements in a same row
	int * restrict col; // Column index of elemetns.
	Carray vec;         // column oriented vector.
};

typedef struct CCS * restrict CCSmat;



/**** Compressed Column storage of an n x n Real Matrix ****/

struct RCCS{
    int  m; // max number of non zero elements in a same row
	int * restrict col; // Column index of elemetns.
	Rarray vec;         // column oriented vector.
};

typedef struct RCCS * restrict RCCSmat;

#endif
