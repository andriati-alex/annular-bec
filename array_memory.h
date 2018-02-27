#ifndef _array_memory_h
#define _array_memory_h

// Default C library
#include <stdio.h>
#include <stdlib.h>

#include "array.h" // complex matrix and vectors datatype

// Memory allocation of complex vector
Rarray rarrDef(unsigned int n);

// Memory allocation of complex vector
Carray carrDef(unsigned int n);

// Memory allocation of real matrix
Rmatrix rmatDef(unsigned int m, unsigned int n);

// Memory allocation of complex matrix
Cmatrix cmatDef(unsigned int m, unsigned int n);

// Release memory of allocated real vector
void rarrFree(Rarray v);

// Release memory of allocated complex vector
void carrFree(Carray v);

// Release memory of allocated complex matrix
void rmatFree(unsigned int m, Rmatrix M);

// Release memory of allocated complex matrix
void cmatFree(unsigned int m, Cmatrix M);

// Printer functions. Vectors as a column matrix
void cPrint(double complex z);
void carrPrint(unsigned int n, Carray v);
void cmatPrint(unsigned int m, unsigned int n, Cmatrix M);

#endif
