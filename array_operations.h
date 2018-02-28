#ifndef _array_operations_h
#define _array_operations_h

// C default library
#include <math.h>

/* Array Data-type */
#include "array.h"

// Fill the elements of an array with a constant number z
void carrFill(unsigned int n, double complex z, Carray v);

// Copy elements "from" array "to" array
void carrCopy(unsigned int n, Carray from, Carray to);

// Real part vector
void carrRPart(unsigned int n, Carray v, Rarray vreal);

// Imaginary part vector
void carrIPart(unsigned int n, Carray v, Rarray vimag);

// Each position addition
void carrAdd(unsigned int n, Carray v1, Carray v2, Carray v);

// convention v1 - v2
void carrSub(unsigned int n, Carray v1, Carray v2, Carray v);

// Each position product
void carrProd(unsigned int n, Carray v1, Carray v2, Carray v);

// multiply all elements by a complex scalar
void carrScalar(unsigned int n, Carray v, double complex z, Carray ans);

// Each position division
void carrDiv(unsigned int n, Carray v1, Carray v2, Carray v);

// Update-like operation v[i] = v1[i] + z * v2[i];
void carrMix(unsigned int n, double complex z, Carray v1, Carray v2, Carray v);

// Each position conjugation
void carrConj(unsigned int n, Carray v, Carray v_conj);

// Take absolute(complex) value for each component
void carrAbs(unsigned int n, Carray v, Rarray vabs);

// Take absolute(complex) squared value for each component
void carrAbs2(unsigned int n, Carray v, Rarray vabs);

// Default scalar product. Convention <v1*, v2>
double complex carrDot(unsigned int n, Carray v1, Carray v2);

// Product and sum elements. Convention <v1, v2>
double complex carrDot2(unsigned int n, Carray v1, Carray v2);

// Vector modulus
double carrMod(unsigned int n, Carray v);

// Vector Modulus squared
double carrMod2(unsigned int n, Carray v);



/************************** ELEMENTARY FUNCTIONS **************************/



// take exponential of each component of z * v.
void carrExp(unsigned int n, double complex z, Carray v, Carray ans);

#endif
