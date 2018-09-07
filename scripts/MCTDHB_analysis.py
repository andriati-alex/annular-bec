
import numpy as np;
import math;
import cmath;
from numba import jit, prange, uint32, int32, float64, complex128;

@jit(uint32(uint32), nopython=True, nogil=True)
def fac(n):
    nfac = 1;
    for i in prange(2, n + 1): nfac = nfac * i;
    return nfac;

@jit(uint32(uint32, uint32), nopython=True, nogil=True)
def NC(Npar, Morb):
    n = 1;
    for i in prange(Npar + Morb - 1, Npar, -1): n = n * i;
    return n / fac(Morb - 1);

@jit((int32, int32, int32, int32[:]), nopython=True, nogil=True)
def IndexToFock(k, N, M, v):
    x = 0;
    m = M - 1;
    for i in prange(0, M, 1): v[i] = 0;
    while (k > 0):
        while (k - NC(N, m) < 0): m = m - 1;
        x = k - NC(N, m);
        while (x >= 0):
            v[m] = v[m] + 1;
            N = N - 1;
            k = x;
            x = x - NC(N, m);
    for i in prange(N, 0, -1): v[0] = v[0] + 1;

@jit((int32, int32, int32[:,:]), nopython=True, nogil=True)
def MountNCmat(N, M, NCmat):
    for i in prange(0, N + 1, 1):
        NCmat[i][0] = 0;
        for j in prange(1, M + 1, 1): NCmat[i][j] = NC(i, j);

@jit((int32, int32, int32[:,:]), nopython=True, nogil=True)
def MountFocks(N, M, IF):
    for k in prange(0, NC(N, M), 1): IndexToFock(k, N, M, IF[k]);

def GetNCmat(N, M):
    NCmat = np.empty([N + 1, M + 1], dtype=np.int32);
    MountNCmat(N, M, NCmat);
    return NCmat;

def GetFocks(N, M):
    IF = np.empty([NC(N, M), M], dtype=np.int32);
    MountFocks(N, M, IF);
    return IF;

@jit(int32(int32, int32, int32[:,:], int32[:]), nopython=True, nogil=True)
def FockToIndex(N, M, NCmat, v):
    n = 0;
    k = 0;
    for i in prange(M - 1, 0, -1):
        n = v[i]; # Number of particles in the orbital
        while (n > 0):
            k = k + NCmat[N][i]; # number of combinations needed
            N = N - 1;           # decrease the number of particles
            n = n - 1;
    return k;

@jit((int32, int32, int32[:,:], int32[:,:], complex128[:], complex128[:,:]),
      nopython=True, nogil=True)

def OBrho(N, M, NCmat, IF, C, rho):
    
    # Initialize variables
    j = 0;
    mod2 = 0.0;
    nc = NCmat[N][M];
    RHO = 0.0 + 0.0j;

    for k in prange(0, M, 1):
        RHO = 0.0 + 0.0j;
        for i in prange(0, nc, 1):
            mod2 = C[i].real * C[i].real + C[i].imag * C[i].imag;
            RHO += mod2 * IF[i][k];
        rho[k][k] = RHO;

        for l in prange(k + 1, M, 1):
            RHO = 0;
            for i in prange(0, nc, 1):
                if (IF[i][k] < 1): continue;
                IF[i][k] -= 1;
                IF[i][l] += 1;
                j = FockToIndex(N, M, NCmat, IF[i]);
                IF[i][k] += 1;
                IF[i][l] -= 1;
                RHO += C[i].conjugate() * C[j] * \
                       math.sqrt((IF[i][l] + 1) * IF[i][k]);
            rho[k][l] = RHO;

    for k in prange(0, M, 1):
        for l in prange(k + 1, M, 1): rho[l][k] = rho[k][l].conjugate();
