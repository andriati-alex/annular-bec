#ifndef _array_h
#define _array_h

#include <complex.h>

/* DATA TYPE FOR LINEAR ALGEBRA MATRICES AND VECTORS
 * *************************************************
 *
 * Pointers declared with restrict keyword  suppose
 * no aliasing between pointers. For this you can't
 * pass the same pointer in two different arguments
 * of the routines.
 *
 * For example don't call carrAdd(n, v, u, v) aiming
 * produce the effect of v[i] = u[i] + v[i].
 *
 * Sometimes it'd be useful to call routines like
 * that way to save memory,  but lose performance
 *
 * *************************************************/




/***************** Real vector and matrix ******************/
typedef double         *  restrict Rarray;
typedef double         ** restrict Rmatrix;

/**************** Complex vector and matrix ****************/
typedef double complex *  restrict Carray;
typedef double complex ** restrict Cmatrix;




/******* Compressed Column storage of an n x n Matrix ******/
struct CCS{
    int  m; // max number of non zero elements in a same row
	int * restrict col; // Column index of elemetns. Size = m * n
	double complex * restrict vec; // column oriented vector. Size = m * n
};

typedef struct CCS * restrict CCSmat;




/**** Compressed Column storage of an n x n Real Matrix ****/
struct RCCS{
    int  m; // mas number of non zero elements in a same row
	int * restrict col; // Column index of elemetns. Size = m * n
	double * restrict vec; // column oriented vector. Size = m * n
};

typedef struct RCCS * restrict RCCSmat;

#endif
