
import math;
import cmath;
import numpy as np;
import scipy.linalg as la;
from numba import jit, prange, uint32, int32, float64, complex128;

"""

    Module devoted to analyze results from MCTDHB method
    ----------------------------------------------------

    The functions contained in this module support the analysis of results 
    from (imaginary/real)time propagation.

"""





@jit(uint32(uint32), nopython=True, nogil=True)

def fac(n):
    """ return n! """
    nfac = 1;
    for i in prange(2, n + 1): nfac = nfac * i;
    return nfac;





@jit(uint32(uint32, uint32), nopython=True, nogil=True)

def NC(Npar, Morb):
    """ return (Npar + Morb - 1)! / ( (Npar)! x (Morb - 1)! )"""
    n = 1;
    for i in prange(Npar + Morb - 1, Npar, -1): n = n * i;
    return n / fac(Morb - 1);





@jit((int32, int32, int32, int32[:]), nopython=True, nogil=True)
def IndexToFock(k, N, M, v):
    """
    Calling: (void) IndexToFock(k, N, M, v)
    --------

    Arguments:
    ----------
    k : Index of occupation configuration coeficients
    N : # of particles
    M : # of orbitals
    v : End up with occupation vector(Fock state) of length Morb
    """
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
    """ Auxiliar of GetNCmat """
    for i in prange(0, N + 1, 1):
        NCmat[i][0] = 0;
        for j in prange(1, M + 1, 1): NCmat[i][j] = NC(i, j);





@jit((int32, int32, int32[:,:]), nopython=True, nogil=True)

def MountFocks(N, M, IF):
    """ Auxiliar of GetFocks """
    for k in prange(0, NC(N, M), 1): IndexToFock(k, N, M, IF[k]);





def GetNCmat(N, M):
    """
    Calling: (numpy 2D array of ints) = GetNCmat(N, M)
    --------
    Returned matrix NC(n,m) = (n + m - 1)! /( (n)! (m - 1)! ).

    arguments:
    ----------
    N: # of particles
    M: # of orbitals
    """
    NCmat = np.empty([N + 1, M + 1], dtype=np.int32);
    MountNCmat(N, M, NCmat);
    return NCmat;





def GetFocks(N, M):
    """
    Calling : (numpy 2D array of ints) = GetFocks(N, M)
    ---------
    Row k has the occupation vector corresponding to C[k].

    arguments:
    ----------
    N : # of particles
    M : # of orbitals
    """
    IF = np.empty([NC(N, M), M], dtype=np.int32);
    MountFocks(N, M, IF);
    return IF;





@jit(int32(int32, int32, int32[:,:], int32[:]), nopython=True, nogil=True)

def FockToIndex(N, M, NCmat, v):
    """
    Calling: (int) = FockToIndex(N, M, NCmat, v)
    --------
    k = Index of Fock Configuration Coeficient of v

    arguments:
    ----------
    N : # of particles
    M : # of orbitals
    NCmat : see GetNCmat function
    """
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
    """
    Calling: (void) OBrho(N, M, NCmat, IF, C, rho)
    --------
    Setup the rho argument with one-body densit matrix

    arguments:
    ----------
    N : # of particles
    M : # of orbitals
    NCmat : see function GetNCmat
    IF : see function GetFocks
    C : coeficients of Fock-configuration states
    rho : Empty M x M matrix. End up configured with values
    """
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





@jit((int32, float64[:], complex128[:,:]), nopython=True, nogil=True)

def EigSort(Nvals, RHOeigvals, RHOeigvecs):
    """
    Calling: (void) EigSort(Nvals, RHOeigvals, RHOeigvecs)
    --------
    Sort the order of eigenvalues to be decreasing
    and the order of columns of eigenvectors.

    arguments
    ---------
    Nvals : dimension of rho = # of orbitals
    RHOeigvals : End up with eigenvalues in decreasing order
    RHOeigvecs : Change the order of columns accordingly to eigenvalues
    """
    auxR = 0.0;
    auxC = 0.0;
    for i in prange(1, Nvals, 1):
        j = i;
        while (RHOeigvals[j] > RHOeigvals[j-1] and j > 0):
            # Sort the vector
            auxR = RHOeigvals[j-1];
            RHOeigvals[j-1] = RHOeigvals[j];
            RHOeigvals[j] = auxR;
            # Sort the matrix
            for k in prange(0, Nvals, 1):
                auxC = RHOeigvecs[k][j-1];
                RHOeigvecs[k][j-1] = RHOeigvecs[k][j];
                RHOeigvecs[k][j] = auxC;
            j = j - 1;





def NatOrb(RHOeigvecs, Orb):
    """
    CALLING:
    --------
    ( 2D numpy array [Morb x Mpos] ) = NatOrb(RHOeigvecs, Orb)

    Arguments:
    ----------
    RHOeigvecs : Matrix with eigenvectors of rho in columns
    Orb : 2D array [Morb x Mpos] given by time propagation
    """
    return np.matmul(RHOeigvecs.conj().T, Orb);




def VonNeumannS(RHOeigvals):
    N = round(RHOeigvals.sum().real);
    return - ( (RHOeigvals.real / N) * np.log(RHOeigvals.real / N) ).sum()
