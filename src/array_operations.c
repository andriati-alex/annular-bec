#include "../include/array_operations.h"



/*          ***********************************************          */
/*                 SETUP VALUES IN VECTOR COMPONENTS                 */
/*          ***********************************************          */



void carrFill(int n, double complex z, Carray v)
{
    int i;
    for (i = 0; i < n; i++) { v[i] = z; }
}

void rarrFill(int n, double x, Rarray v)
{
    int i;
    for (i = 0; i < n; i++) { v[i] = x; }
}

void rarrFillInc(int n, double x0, double dx, Rarray v)
{   // Fill increasing along the elements with a step dx
    // starting from initial value x0
    int i;
    v[0] = x0;
    for (i = 1; i < n; i++) v[i] = v[i - 1] + dx;
}

void carrCopy(int n, Carray from, Carray to)
{ int i; for (i = 0; i < n; i++) to[i] = from[i]; }

void rarrCopy(int n, Rarray from, Rarray to)
{ int i; for (i = 0; i < n; i++) to[i] = from[i]; }

void fromMKL(int n, CMKLarray a, Carray b)
{   // Copy data from MKL array to Complex array
    int i;
    for (i = 0; i < n; i++) b[i] = a[i].real + I * a[i].imag;
}

void toMKL(int n, CMKLarray a, Carray b)
{   // Copy data from Complex array to MKL datatype
    int i;
    for (i = 0; i < n; i++)
    {
        a[i].real = creal(b[i]);
        a[i].imag = cimag(b[i]);
    }
}



/*          ***********************************************          */
/*                   BASIC ELEMENT-WISE OPERATIONS                   */
/*          ***********************************************          */



void carrRPart(int n, Carray v, Rarray vreal)
{ int i; for (i = 0; i < n; i++) vreal[i] = creal(v[i]); }

void carrIPart(int n, Carray v, Rarray vimag)
{ int i; for (i = 0; i < n; i++) vimag[i] = cimag(v[i]); }

void carrConj(int n, Carray v, Carray v_conj)
{ int i; for (i = 0; i < n; i++) v_conj[i] = conj(v[i]); }

void carrAdd(int n, Carray v1, Carray v2, Carray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] + v2[i]; }

void rarrAdd(int n, Rarray v1, Rarray v2, Rarray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] + v2[i]; }

void carrSub(int n, Carray v1, Carray v2, Carray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] - v2[i]; }

void rarrSub(int n, Rarray v1, Rarray v2, Rarray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] - v2[i]; }

void carrMultiply(int n, Carray v1, Carray v2, Carray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] * v2[i]; }

void rarrMultiply(int n, Rarray v1, Rarray v2, Rarray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] * v2[i]; }

void carrScalarMultiply(int n, Carray v, double complex z, Carray ans)
{ int i; for (i = 0; i < n; i++) ans[i] = v[i] * z; }

void rarrScalarMultiply(int n, Rarray v, double z, Rarray ans)
{ int i; for (i = 0; i < n; i++) ans[i] = v[i] * z; }

void carrScalarAdd(int n, Carray v, double complex z, Carray ans)
{ int i; for (i = 0; i < n; i++) ans[i] = v[i] + z; }

void rarrScalarAdd(int n, Rarray v, double z, Rarray ans)
{ int i; for (i = 0; i < n; i++) ans[i] = v[i] + z; }

void carrDiv(int n, Carray v1, Carray v2, Carray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] / v2[i]; }

void rarrDiv(int n, Rarray v1, Rarray v2, Rarray v)
{ int i; for (i = 0; i < n; i++) v[i] = v1[i] / v2[i]; }

void carrUpdate(int n, Carray v1, double complex z, Carray v2, Carray v)
{
    int i;
    for (i = 0; i < n; i++) v[i] = v1[i] + z * v2[i];
}

void rcarrUpdate(int n, Carray v1, double complex z, Rarray v2, Carray v)
{
    int i;
    for (i = 0; i < n; i++) v[i] = v1[i] + z * v2[i];
}

void rarrUpdate(int n, Rarray v1, double z, Rarray v2, Rarray v)
{
    int i;
    for (i = 0; i < n; i++) v[i] = v1[i] + z * v2[i];
}

void carrAbs(int n, Carray v, Rarray vabs)
{ int i; for (i = 0; i < n; i++) vabs[i] = cabs(v[i]); }

void rarrAbs(int n, Rarray v, Rarray vabs)
{ int i; for (i = 0; i < n; i++) vabs[i] = fabs(v[i]); }

void rarrAbs2(int n, Rarray v, Rarray vabs)
{ int i; for (i = 0; i < n; i++) vabs[i] = v[i] * v[i]; }

void carrAbs2(int n, Carray v, Rarray vabs)
{
    int i;
    for (i = 0; i < n; i++)
    {
        vabs[i] = creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }
}

void renormalizeVector(int n, Carray v, double norm)
{
    int i;
    double renorm;

    renorm = norm / carrMod(n, v);
    for (i = 0; i < n; i ++) v[i] = v[i] * renorm;
}



/*          ***********************************************          */
/*                             FUNCTIONALS                           */
/*          ***********************************************          */



double complex carrDot(int n, Carray v1, Carray v2)
{
    double complex z = 0;
    int i;
    for (i = 0; i < n; i++) z += conj(v1[i]) * v2[i];
    return z;
}

double complex carrDot2(int n, Carray v1, Carray v2)
{
    double complex z = 0;
    int i;
    for (i = 0; i < n; i++) { z += v1[i] * v2[i]; }
    return z;
}

double rarrDot(int n, Rarray v1, Rarray v2)
{
    double z = 0;
    int i;
    for (i = 0; i < n; i++) { z += v1[i] * v2[i]; }
    return z;
}

double carrMod(int n, Carray v)
{
    double mod = 0;
    int i;
    for (i = 0; i < n; i++)
    {
        mod += creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }
    return sqrt(mod);
}

double carrMod2(int n, Carray v)
{
    double mod = 0;
    int i;
    for (i = 0; i < n; i++)
    {
        mod += creal(v[i]) * creal(v[i]) + cimag(v[i]) * cimag(v[i]);
    }
    return mod;
}

double complex carrReduction(int n, Carray v)
{
    double complex red = 0;
    int i;
    for (i = 0; i < n; i++) red += v[i];
    return red;;
}

double rarrReduction(int n, Rarray v)
{
    double red = 0;
    int i;
    for (i = 0; i < n; i++) red += v[i];
    return red;;
}



/*          ***********************************************          */
/*                   FUNCTION COMPUTED ELEMENT-WISE                  */
/*          ***********************************************          */



void carrExp(int n, double complex z, Carray v, Carray ans)
{
    int i;
    #pragma omp parallel for private(i)
    for (i = 0; i < n; i++) { ans[i] = cexp(z * v[i]); }
}

void rcarrExp(int n, double complex z, Rarray v, Carray ans)
{
    int i;
    #pragma omp parallel for private(i)
    for (i = 0; i < n; i++) { ans[i] = cexp(z * v[i]); }
}
