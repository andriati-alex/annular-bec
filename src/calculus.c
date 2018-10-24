#include "../include/calculus.h"





double complex Csimps(int n, Carray f, double dx)
{
    if (n == 1) return 0; // integrating over empty set

    int i, m = n - 1;
    double complex sums = 0;

    for (i = 0; i < m / 2 - 1; i++)
    {
        sums = sums + 2 * f[2*(i+1)] + 4 * f[2*i+1];
    }

    if (m % 2 == 0) { return (sums + f[0] + f[m] + 4 * f[m-1]) * dx / 3; }

    // case odd number of intervals integrate last by trapezium
    if (m > 2) sums = (sums + f[0] + f[m-1] + 4 * f[m-2]) * dx / 3;
    return sums + (f[m-1] + f[m]) * dx * 0.5;
}





double Rsimps(int n, Rarray f, double dx)
{
    if (n == 1) return 0; // integrating over empty set

    int i, m = n - 1;
    double sums = 0;

    for (i = 0; i < m / 2 - 1; i++) sums += 2 * f[2*(i+1)] + 4 * f[2*i+1];

    if (m % 2 == 0) return (sums + f[0] + f[m] + 4 * f[m-1]) * dx / 3;

    /* case odd number of intervals integrate last by trapezium */
    if (m > 2) sums = (sums + f[0] + f[m-1] + 4 * f[m-2]) * dx / 3;
    return sums + (f[m-1] + f[m]) * dx * 0.5;
}





void renormalize(int n, Carray f, double dx, double norm)
{
    int i;
    double renorm;
    Rarray ToInt = rarrDef(n);

    carrAbs2(n, f, ToInt);
    renorm = norm * sqrt(1.0 / Rsimps(n, ToInt, dx));
    for (i = 0; i < n; i++) f[i] = f[i] * renorm;

    free(ToInt);
}





void Ortonormalize(int Mfun, int Mpos, double dx, Cmatrix F)
{   // F[k][:] has the k-th function of the basis
    int i, j, k;
    Carray toInt = carrDef(Mpos);

    renormalize(Mpos, F[0], dx, 1.0);
    for (i = 1; i < Mfun; i++)
    {
        for (j = 0; j < i; j++)
        {   // The projection are integrals of the product below
            for (k = 0; k < Mpos; k++)
                toInt[k] = conj(F[j][k]) * F[i][k];
            // Iterative Gram-Schmidt (see wikipedia)
            for (k = 0; k < Mpos; k++)
                F[i][k] = F[i][k] - Csimps(Mpos, toInt, dx) * F[j][k];
        }
        // normalized to unit the new vector
        renormalize(Mpos, F[i], dx, 1.0);
    }

    free(toInt);
}





void dxFFT(int n, Carray f, double dx, Carray dfdx)
{
    int i,
        N;

    N = n - 1; // Assumes the connection f[n-1] = f[0] the boundary

    double freq;

    carrCopy(N, f, dfdx);

    MKL_LONG s; // status of called MKL FFT functions

    DFTI_DESCRIPTOR_HANDLE desc;
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, N);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt((double) N));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt((double) N));
    // s = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    s = DftiCommitDescriptor(desc);

    s = DftiComputeForward(desc, dfdx);

    for (i = 0; i < N; i++) {
        if (i <= (N - 1) / 2) { freq = (2 * PI * i) / (N * dx);       }
        else                  { freq = (2 * PI * (i - N)) / (N * dx); }
        dfdx[i] *= freq * I;
    }

    s = DftiComputeBackward(desc, dfdx);
    s = DftiFreeDescriptor(&desc);
    dfdx[N] = dfdx[0]; // boundary point
}





void dxCyclic(int n, Carray f, double dx, Carray dfdx)
{
    int i;
    double r = 1.0 / (2 * dx);

    dfdx[0]   = (f[1] - f[n-2]) * r;
    dfdx[n-2] = (f[0] - f[n-3]) * r;
    dfdx[n-1] = dfdx[0];
    for (i = 1; i < n - 1; i++) dfdx[i] = (f[i+1] - f[i-1]) * r;
}





void applyL0(int n, Carray f, double dx, double a2, double complex a1, 
             Rarray V, double inter, double mu, Carray L0f)
{
    int i;
    double complex ddx = a1 / (2 * dx);
    double abs2, d2dx = a2 / (dx * dx);

    #pragma omp parallel for private(i, abs2)
    for (i = 1; i < n - 1; i++) {
        abs2 = creal(f[i]) * creal(f[i]) + cimag(f[i]) * cimag(f[i]);
        L0f[i] =  (f[i + 1] - f[i - 1]) * ddx;
        L0f[i] += (f[i + 1] - 2 * f[i] + f[i - 1]) * d2dx;
        L0f[i] += f[i] * (V[i] - mu + inter * abs2);
    }
    
    abs2 = creal(f[0]) * creal(f[0]) + cimag(f[0]) * cimag(f[0]);
    L0f[0]  =  (f[1] - f[n - 1]) * ddx;
    L0f[0] += (f[1] - 2 * f[0] + f[n - 1]) * d2dx;
    L0f[0] += f[0] * (V[0] - mu + inter * abs2);
    
    abs2 = creal(f[n-1]) * creal(f[n-1]) + cimag(f[n-1]) * cimag(f[n-1]);
    L0f[n-1]  =  (f[0] - f[n - 2]) * ddx;
    L0f[n-1] += (f[0] - 2 * f[n-1] + f[n - 2]) * d2dx;
    L0f[n-1] += f[n-1] * (V[n-1] - mu + inter * abs2);
}





double complex Functional(int M, double dx, double a2, double complex a1,
               double inter, Rarray V, Carray f)
{
    int i;
    double norm;
    double complex E;

    Carray DF = carrDef(M);
    Rarray abs2F = rarrDef(M);
    Rarray abs2DF = rarrDef(M);
    Carray Int = carrDef(M);

    #pragma omp parallel sections firstprivate(dx, inter)
    {
        #pragma omp section
        {
            dxCyclic(M, f, dx, DF);
            carrAbs2(M, DF, abs2DF);
        }
        #pragma omp section
        {
            carrAbs2(M, f, abs2F);
            for (i = 0; i < M; i++)
                Int[i] = (V[i] + inter * abs2F[i]) * abs2F[i];
        }
    }

    #pragma parallel for private(i)
    for (i = 0; i < M; i++)
        Int[i] += - a2 * abs2DF[i] + a1 * conj(f[i]) * DF[i];

    #pragma omp parallel sections
    {
        #pragma omp section
        { E = Csimps(M, Int, dx);      }
        #pragma omp section
        { norm = Rsimps(M, abs2F, dx); }
    }

    // release memory
    free(DF); free(abs2F); free(abs2DF); free(Int);

    return E / norm;
}
