#ifndef _array_memory_h
#define _array_memory_h

#include <stdio.h>
#include <stdlib.h>
#include "array.h"



Rarray rarrDef(int n);
// Allocate real vector

Carray carrDef(int n);
// Allocate complex vector

CMKLarray CMKLdef(int n);
// Allocate MKL's complex vector

Rmatrix rmatDef(int m, int n);
// Allocate real matrix with m rows and n columns

Cmatrix cmatDef(int m, int n);
// Allocate complex matrix with m rows and n columns

CCSmat ccsmatDef(int n, int max_nonzeros);
// Allocate CCS matrix structure with n rows



void rmatFree(int m, Rmatrix M);
// Release allocated real matrix with m rows

void cmatFree(int m, Cmatrix M);
// Release allocated complex matrix with m rows

void CCSFree(CCSmat M);
// Relase CCS matrix

void RCCSFree(RCCSmat M);
// Relase CCS matrix with real entries

#endif
