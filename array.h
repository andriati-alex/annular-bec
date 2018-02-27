#ifndef _array_h
#define _array_h

#include <complex.h>

typedef double         *Rarray;
typedef Rarray         *Rmatrix;
typedef double complex *Carray;
typedef Carray         *Cmatrix;

struct CCS{   // Compressed Column storage of an n x n Matrix
	int  m;   // mas number of non zero elements in a same row
	int *col; // Column index of elemetns. Size = m * n
	double complex *vec; // column oriented vector. Size = m * n
};  typedef struct CCS * CCSmat;

#endif
