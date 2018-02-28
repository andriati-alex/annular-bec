/****** Header file ******/
#include "array_operations.h"

void carrFill(unsigned int n, double complex z, Carray v)
{ for (unsigned int i = 0; i < n; i++) { v[i] = z; } }

void carrCopy(unsigned int n, Carray from, Carray to)
{ for (unsigned int i = 0; i < n; i++) { to[i] = from[i]; } }

void carrRPart(unsigned int n, Carray v, Rarray vreal)
{ for (unsigned int i = 0; i < n; i++) { vreal[i] = creal(v[i]); } }

void carrIPart(unsigned int n, Carray v, Rarray vimag)
{ for (unsigned int i = 0; i < n; i++) { vimag[i] = cimag(v[i]); } }

void carrAdd(unsigned int n, Carray v1, Carray v2, Carray v)
{ for (unsigned int i = 0; i < n; i++) { v[i] = v1[i] + v2[i]; } }

void carrSub(unsigned int n, Carray v1, Carray v2, Carray v)
{ for (unsigned int i = 0; i < n; i++) { v[i] = v1[i] - v2[i]; } }

void carrProd(unsigned int n, Carray v1, Carray v2, Carray v)
{ for (unsigned int i = 0; i < n; i++) { v[i] = v1[i] * v2[i]; } }

void carrScalar(unsigned int n, Carray v, double complex z, Carray ans)
{ for (unsigned int i = 0; i < n; i++) { ans[i] = v[i] * z; } }

void carrDiv(unsigned int n, Carray v1, Carray v2, Carray v)
{ for (unsigned int i = 0; i < n; i++) { v[i] = v1[i] / v2[i]; } }

void carrMix(unsigned int n, double complex z, Carray v1, Carray v2, Carray v)
{ for (unsigned int i = 0; i < n; i++) { v[i] = v1[i] + z * v2[i]; } }

void carrConj(unsigned int n, Carray v, Carray v_conj)
{ for (unsigned int i = 0; i < n; i++) { v_conj[i] = conj(v[i]); } }

void carrAbs(unsigned int n, Carray v, Rarray vabs)
{ for (unsigned int i = 0; i < n; i++) { vabs[i] = cabs(v[i]); } }

void carrAbs2(unsigned int n, Carray v, Rarray vabs)
{
    for (unsigned int i = 0; i < n; i++) {
        vabs[i] = creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }
}

double complex carrDot(unsigned int n, Carray v1, Carray v2) {
    double complex z = 0 + 0 * I;
    for (unsigned int i = 0; i < n; i++) { z += conj(v1[i]) * v2[i]; }
    return z;
}

double complex carrDot2(unsigned int n, Carray v1, Carray v2) {
    unsigned int i;
    double complex z = 0 + 0 * I;
    for (i = 0; i < n; i++) { z += v1[i] * v2[i]; }
    return z;
}

double carrMod(unsigned int n, Carray v) {
    double mod = 0;
    unsigned int i;
    for (i = 0; i < n; i++) {
        mod += creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }
    return sqrt(mod);
}

double carrMod2(unsigned int n, Carray v) {
    double mod = 0;
    unsigned int i;
    for (i = 0; i < n; i++) {
        mod += creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }
    return mod;
}



/************************** ELEMENTARY FUNCTIONS **************************/



void carrExp(unsigned int n, double complex z, Carray v, Carray ans)
{ for (unsigned int i = 0; i < n; i++) { ans[i] = cexp(z * v[i]); } }
