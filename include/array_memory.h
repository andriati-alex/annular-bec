#ifndef _array_memory_h
#define _array_memory_h

#include <stdio.h>
#include <stdlib.h>
#include "array.h"



/* ======================================================================== */
/* ======================== ALLOCATE MEMORY SECTION ======================= */
/* ======================================================================== */



// Allocate real vector
Rarray rarrDef(int n);

// Allocate complex vector
Carray carrDef(int n);

// Allocate MKL's complex vector
CMKLarray CMKLdef(int n);

// Allocate real matrix with m rows and n columns
Rmatrix rmatDef(int m, int n);

// Allocate complex matrix with m rows and n columns
Cmatrix cmatDef(int m, int n);



/* ======================================================================== */
/* ======================== RELEASE MEMORY SECTION ======================== */
/* ======================================================================== */



// Release allocated real vector
void rarrFree(Rarray v);

// Release allocated complex vector
void carrFree(Carray v);

// Release allocated real matrix with m rows
void rmatFree(int m, Rmatrix M);

// Release allocated complex matrix with m rows
void cmatFree(int m, Cmatrix M);

// Relase CCS matrix
void CCSFree(CCSmat M);
void RCCSFree(RCCSmat M);



/* ======================================================================== */
/* ============================  PRINT SECTION ============================ */
/* ======================================================================== */



void cPrint(double complex z);

void carrPrint(int n, Carray v); // as column matrix

void rarrPrint(int n, Rarray v); // as column matrix

void cmat_print(int m, int n, Cmatrix M);

// Print array in a text file suitable to import with python
void carr_txt(char fname [], int M, Carray v);
    
// Print matrix in a text file suitable to import with python
void cmat_txt(char fname [],
              int N, int row_step, int M, int col_step, Cmatrix S);

#endif
