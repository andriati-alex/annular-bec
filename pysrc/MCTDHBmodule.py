
import cmath;
import numpy as np;
import scipy.linalg as la;

from math import sqrt, pi;
from scipy.integrate import simps;
from numba import jit, prange, int32, uint32, uint64, int64, float64, complex128;





"""
=============================================================================


    Module devoted to analyze results from MCTDHB method
    ----------------------------------------------------


    The functions contained in this module support the analysis of results 
    from (imaginary/real)time propagation. Some of the functions  are  the
    same as those used in C language, and use NUMBA compilation to improve
    performance


=============================================================================
"""





@jit([ uint32(uint32) , uint64(uint64) ], nopython=True, nogil=True)

def fac(n):
    """ return n! """
    nfac = 1;
    for i in prange(2, n + 1): nfac = nfac * i;
    return nfac;





@jit( uint32(uint32, uint32), nopython=True, nogil=True)

def NC(Npar, Morb):
    """ return (Npar + Morb - 1)! / ( (Npar)! x (Morb - 1)! )"""
    n = 1;
    if (Npar > Morb):
        for i in prange(Npar + Morb - 1, Npar, -1): n = n * i;
        return n / fac(Morb - 1);
    else :
        for i in prange(Npar + Morb - 1, Morb - 1, -1): n = n * i;
        return n / fac(Npar);





@jit( [ (int32, int32, int32, int32[:]), (int64, int32, int32, int32[:]) ],
      nopython=True, nogil=True)

def IndexToFock(k, N, M, v):
    """
    Calling: (void) IndexToFock(k, N, M, v)
    -------

    Arguments:
    ---------
    k : Index of configuration-state coefficient
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





@jit( int32(int32, int32, int32[:,:], int32[:]) , nopython=True, nogil=True)

def FockToIndex(N, M, NCmat, v):
    """
    Calling: (int) = FockToIndex(N, M, NCmat, v)
    --------
    k = Index of Fock Configuration Coeficient of v

    arguments:
    ----------
    N     : # of particles
    M     : # of orbitals
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





@jit( [ (int32, int32, int32[:,:]), (int32, int32, int64[:,:]) ],
      nopython=True, nogil=True)

def MountNCmat(N, M, NCmat):
    """ Auxiliar of GetNCmat """
    for i in prange(0, N + 1, 1):
        NCmat[i][0] = 0;
        for j in prange(1, M + 1, 1): NCmat[i][j] = NC(i, j);





@jit( (int32, int32, int32[:,:]) , nopython=True, nogil=True)

def MountFocks(N, M, IF):
    """ Auxiliar of GetFocks """
    for k in prange(0, NC(N, M), 1): IndexToFock(k, N, M, IF[k]);





def GetNCmat(N, M):
    """
    Calling: (numpy 2D array of ints) = GetNCmat(N, M)
    --------
    Returned matrix NC(n,m) = (n + m - 1)! / ( n! (m - 1)! ).

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
    -------
    Row k has the occupation vector corresponding to C[k].

    arguments :
    ---------
    N : # of particles
    M : # of orbitals
    """
    IF = np.empty([NC(N, M), M], dtype=np.int32);
    MountFocks(N, M, IF);
    return IF;





@jit( (int32, int32, int32[:,:], int32[:,:], complex128[:], complex128[:,:]),
      nopython=True, nogil=True)

def OBrho(N, M, NCmat, IF, C, rho):

    """
    Calling : (void) OBrho(N, M, NCmat, IF, C, rho)
    -------
    Setup rho argument with one-body densit matrix

    arguments :
    ---------
    N     : # of particles
    M     : # of orbitals (also dimension of rho)
    NCmat : see function GetNCmat
    IF    : see function GetFocks
    C     : coeficients of Fock-configuration states
    rho   : Empty M x M matrix. End up configured with values
    """

    # Initialize variables
    j    = 0;
    mod2 = 0.0;
    nc   = NCmat[N][M];
    RHO  = 0.0 + 0.0j;

    for k in prange(0, M, 1):

        RHO = 0.0 + 0.0j;

        for i in prange(0, nc, 1):
            mod2 = C[i].real * C[i].real + C[i].imag * C[i].imag;
            RHO  = RHO + mod2 * IF[i][k];

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
                RHO = RHO + C[i].conjugate() * C[j] * \
                      sqrt((IF[i][l] + 1) * IF[i][k]);

            rho[k][l] = RHO;
            rho[l][k] = RHO.conjugate();

# End OBrho routine ------------------------------------------------------










@jit( (int32, int32, int32[:,:], int32[:,:], complex128[:], complex128[:]),
      nopython=True, nogil=True)

def TBrho(N, M, NCmat, IF, C, rho):

    """
    Calling : (void) OBrho(N, M, NCmat, IF, C, rho)
    -------
    Setup rho argument with two-body density matrix

    arguments :
    ---------
    N     : # of particles
    M     : # of orbitals (also dimension of rho)
    NCmat : see function GetNCmat
    IF    : see function GetFocks
    C     : coeficients of Fock-configuration states
    rho   : Empty M x M matrix. End up configured with values
    """

    # index result of conversion FockToIndex
    j = int(0);

    # occupation vector
    v = np.empty(M, dtype=np.int32);
    
    mod2 = float(0);   # |Cj| ^ 2
    sqrtOf = float(0); # Factors from the action of creation/annihilation

    RHO = float(0);

    # Auxiliar to memory access of two-body matrix
    M2 = M * M;
    M3 = M * M * M;

    nc = NCmat[N][M];





    # ---------------------------------------------
    # Rule 1: Creation on k k / Annihilation on k k
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        RHO = 0;

        for i in prange(0, nc, 1):
            mod2 = C[i].real * C[i].real + C[i].imag * C[i].imag;
            RHO  = RHO + mod2 * IF[i][k] * (IF[i][k] - 1);

        rho[k + M * k + M2 * k + M3 * k] = RHO;
    # ---------------------------------------------------------------------



    # ---------------------------------------------
    # Rule 2: Creation on k s / Annihilation on k s
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for s in prange(k + 1, M, 1):

            RHO = 0;

            for i in prange(0, nc, 1):
                mod2 = C[i].real * C[i].real + C[i].imag * C[i].imag;
                RHO += mod2 * IF[i][k] * IF[i][s];

            # commutation of bosonic operators is used
            # to fill elements by exchange  of indexes
            rho[k + s * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + k * M3] = RHO;
            rho[k + s * M + s * M2 + k * M3] = RHO;
    # ---------------------------------------------------------------------



    # ---------------------------------------------
    # Rule 3: Creation on k k / Annihilation on q q
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for q in prange(k + 1, M, 1):

            RHO = 0;

            for i in prange(0, nc, 1):
                
                if (IF[i][k] < 2) : continue;

                for t in prange(0, M, 1) : v[t] = IF[i][t];

                sqrtOf = sqrt((v[k] - 1) * v[k] * (v[q] + 1) * (v[q] + 2));
                v[k] -= 2;
                v[q] += 2;
                j = FockToIndex(N, M, NCmat, v);
                RHO += C[i].conjugate() * C[j] * sqrtOf;

            # Use 2-index-'hermiticity'
            rho[k + k * M + q * M2 + q * M3] = RHO;
            rho[q + q * M + k * M2 + k * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # ---------------------------------------------
    # Rule 4: Creation on k k / Annihilation on k l
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for l in prange(k + 1, M, 1):

            RHO = 0;

            for i in prange(0, nc, 1):

                if (IF[i][k] < 2) : continue;

                for t in prange(0, M, 1): v[t] = IF[i][t];

                sqrtOf = (v[k] - 1) * sqrt(v[k] * (v[l] + 1));
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);
                RHO += C[i].conjugate()* C[j] * sqrtOf;

            rho[k + k * M + k * M2 + l * M3] = RHO;
            rho[k + k * M + l * M2 + k * M3] = RHO;
            rho[l + k * M + k * M2 + k * M3] = RHO.conjugate();
            rho[k + l * M + k * M2 + k * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # ---------------------------------------------
    # Rule 5: Creation on k s / Annihilation on s s
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for s in prange(k + 1, M, 1):

            RHO = 0;

            for i in prange(0, nc, 1):

                if (IF[i][k] < 1 or IF[i][s] < 1) : continue;

                for t in prange(0, M, 1): v[t] = IF[i][t];

                sqrtOf = v[s] * sqrt(v[k] * (v[s] + 1));
                v[k] -= 1;
                v[s] += 1;
                j = FockToIndex(N, M, NCmat, v);
                RHO += C[i].conjugate() * C[j] * sqrtOf;

            rho[k + s * M + s * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + s * M3] = RHO;
            rho[s + s * M + s * M2 + k * M3] = RHO.conjugate();
            rho[s + s * M + k * M2 + s * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # -----------------------------------------------------------
    # Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for q in prange(k + 1, M, 1):

            for l in prange(q + 1, M, 1):

                RHO = 0;

                for i in prange(0, nc, 1):
                    
                    if (IF[i][k] < 2) : continue;

                    for t in prange(0, M, 1): v[t] = IF[i][t];

                    sqrtOf = sqrt( v[k] * (v[k] - 1) * 
                             (v[q] + 1) * (v[l] + 1) );
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += C[i].conjugate() * C[j] * sqrtOf;

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = RHO.conjugate();
                rho[q + l * M + k * M2 + k * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # -----------------------------------------------------------
    # Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
    # ---------------------------------------------------------------------
    for q in prange(0, M, 1):
        
        for k in prange(q + 1, M, 1):
            
            for l in prange(k + 1, M, 1):

                RHO = 0;

                for i in prange(0, nc, 1):
                    
                    if (IF[i][k] < 2) : continue;

                    for t in prange(0, M, 1): v[t] = IF[i][t];

                    sqrtOf = sqrt( v[k] * (v[k] - 1) * 
                             (v[q] + 1) * (v[l] + 1) );
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += C[i].conjugate() * C[j] * sqrtOf;

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = RHO.conjugate();
                rho[q + l * M + k * M2 + k * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # -----------------------------------------------------------
    # Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
    # ---------------------------------------------------------------------
    for q in prange(0, M, 1):
        
        for l in prange(q + 1, M, 1):
            
            for k in prange(l + 1, M, 1):

                RHO = 0;

                for i in prange(0, nc, 1):
                    
                    if (IF[i][k] < 2) : continue;

                    for t in prange(0, M, 1): v[t] = IF[i][t];

                    sqrtOf = sqrt( v[k] * (v[k] - 1) * 
                             (v[q] + 1) * (v[l] + 1) );
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += C[i].conjugate() * C[j] * sqrtOf;

                rho[k + k * M + q * M2 + l * M3] = RHO;
                rho[k + k * M + l * M2 + q * M3] = RHO;
                rho[l + q * M + k * M2 + k * M3] = RHO.conjugate();
                rho[q + l * M + k * M2 + k * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # -----------------------------------------------------------
    # Rule 7.0: Creation on k s / Annihilation on s l (s < k < l)
    # ---------------------------------------------------------------------
    for s in prange(0, M, 1):

        for k in prange(s + 1, M, 1):

            for l in prange(k + 1, M, 1):

                RHO = 0;

                for i in prange(0, nc, 1):
                    
                    if (IF[i][k] < 1 or IF[i][s] < 1) : continue;

                    for t in prange(0, M, 1) : v[t] = IF[i][t];
                    
                    sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += C[i].conjugate() * C[j] * sqrtOf;

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = RHO.conjugate();
                rho[s + l * M + s * M2 + k * M3] = RHO.conjugate();
                rho[s + l * M + k * M2 + s * M3] = RHO.conjugate();
                rho[l + s * M + k * M2 + s * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # -----------------------------------------------------------
    # Rule 7.1: Creation on k s / Annihilation on s l (k < s < l)
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for s in prange(k + 1, M, 1):

            for l in prange(s + 1, M, 1):

                RHO = 0;

                for i in prange(0, nc, 1):

                    if (IF[i][k] < 1 or IF[i][s] < 1) : continue;

                    for t in prange(0, M, 1) : v[t] = IF[i][t];

                    sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += C[i].conjugate() * C[j] * sqrtOf;

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = RHO.conjugate();
                rho[s + l * M + s * M2 + k * M3] = RHO.conjugate();
                rho[s + l * M + k * M2 + s * M3] = RHO.conjugate();
                rho[l + s * M + k * M2 + s * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # -----------------------------------------------------------
    # Rule 7.2: Creation on k s / Annihilation on s l (k < l < s)
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for l in prange(k + 1, M, 1):

            for s in prange(l + 1, M, 1):

                RHO = 0;

                for i in prange(0, nc, 1):

                    if (IF[i][k] < 1 or IF[i][s] < 1) : continue;

                    for t in prange(0, M, 1) : v[t] = IF[i][t];

                    sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += C[i].conjugate() * C[j] * sqrtOf;

                rho[k + s * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + s * M2 + l * M3] = RHO;
                rho[s + k * M + l * M2 + s * M3] = RHO;
                rho[k + s * M + l * M2 + s * M3] = RHO;

                rho[l + s * M + s * M2 + k * M3] = RHO.conjugate();
                rho[s + l * M + s * M2 + k * M3] = RHO.conjugate();
                rho[s + l * M + k * M2 + s * M3] = RHO.conjugate();
                rho[l + s * M + k * M2 + s * M3] = RHO.conjugate();
    # ---------------------------------------------------------------------



    # ---------------------------------------------
    # Rule 8: Creation on k s / Annihilation on q l
    # ---------------------------------------------------------------------
    for k in prange(0, M, 1):

        for s in prange(0, M, 1):

            if (s == k) : continue;

            for q in prange(0, M, 1):

                if (q == s or q == k) : continue;

                for l in prange(0, M, 1):

                    RHO = 0;

                    if (l == k or l == s or l == q) : continue;

                    for i in prange(0, nc, 1):
                        
                        if (IF[i][k] < 1 or IF[i][s] < 1) : continue;

                        for t in prange(0, M, 1) : v[t] = IF[i][t];

                        sqrtOf = sqrt(v[k] * v[s] * (v[q] + 1) * (v[l] + 1));
                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += C[i].conjugate() * C[j] * sqrtOf;

                    rho[k + s * M + q * M2 + l * M3] = RHO;

                # Finish l loop
            # Finish q loop
        # Finish s loop
    # Finish k loop
    # ---------------------------------------------------------------------

# END OF ROUTINE  ---------------------------------------------------------










def GetOBrho(Npar, Morb, C):
    """
    CALLING : ( 2D array [Morb, Morb] ) = GetOBrho(Npar, Morb, C)
    -------
    return the one-body density matrix rho[k,l] = < a*_k a_l >

    arguments :
    ---------
    Npar : # of particles
    Morb : # of orbitals (dimension of the matrix returned)
    C    : coefficients of configuration states basis
    """
    rho = np.empty([ Morb , Morb ], dtype=np.complex128);
    NCmat = GetNCmat(Npar, Morb);
    IF = GetFocks(Npar, Morb);
    OBrho(Npar, Morb, NCmat, IF, C, rho);
    return rho;





def GetTBrho(Npar, Morb, C):
    """
    CALLING : ( 1D array [ Morb ^4 ] ) = GetOBrho(Npar, Morb, C)
    -------
    return rho[k + s*Morb + q*Morb^2 + l*Morb^3] = < a*_k a*_s a_q a_l >

    arguments :
    ---------
    Npar : # of particles
    Morb : # of orbitals (dimension of the matrix returned)
    C    : coefficients of configuration states basis
    """
    rho = np.empty(Morb * Morb * Morb * Morb, dtype=np.complex128);
    NCmat = GetNCmat(Npar, Morb);
    IF = GetFocks(Npar, Morb);
    TBrho(Npar, Morb, NCmat, IF, C, rho);
    return rho;





@jit( (int32, float64[:], complex128[:,:]) , nopython=True, nogil=True)

def EigSort(Nvals, RHOeigvals, RHOeigvecs):
    """
    Calling : (void) EigSort(Nvals, RHOeigvals, RHOeigvecs)
    -------
    Sort the order of eigenvalues to be decreasing and the order
    of columns of eigenvectors accordingly so that the  k column
    keep being the eigenvector of k-th eigenvalue.

    arguments :
    ---------
    Nvals      : dimension of rho = # of orbitals
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
    CALLING :
    -------
    ( 2D numpy array [Morb x Mpos] ) = NatOrb(RHOeigvecs, Orb)

    Arguments :
    ---------
    RHOeigvecs : Matrix with eigenvectors of rho in columns
    Orb : 2D array [Morb x Mpos] given by time propagation
    """
    return np.matmul(RHOeigvecs.conj().T, Orb);





def TimeOccupation(Morb, Nsteps, rhotime):
    """
    CALLING :
    -------
    ( 2D numpy array [Nsteps, Morb] ) = TimeOccupation(Morb, Nsteps, rhotime)

    arguments :
    ---------
    Morb    : # of orbitals (# of columns in rhotime)
    Nsteps  : # number of time steps (# of rows in rhotime)
    rhotime : each line has row-major matrix representation (a vector)
              that is the  one-body  densiity matrix at each time-step
    """
    eigval = np.empty( [Nsteps, Morb] , dtype=np.complex128 );
    for i in range(Nsteps):
        eigval[i], eigvec = la.eig( rhotime[i].reshape(Morb, Morb) );
        EigSort(Morb, eigval[i].real, eigvec);
    return eigval;





def SpatialOBdensity(M, NOoccu, NO):
    """
    CALLING :
    -------
    ( 2d array [M,M] ) = SpatialOBdensity(M, occu, NO)

    arguments :
    ---------
    M      : # of discrete positions
    NOoccu : occu[k] has # of particles occupying NO[k,:] orbital
    NO     : [Morb,M] matrix with each row being a natural orbital
    """
    n = np.zeros([M , M], dtype=np.complex128);
    for i in range(M):
        for j in range(M):
            for k in range(occu.size):
                n[i,j] = n[i,j] + NOoccu[k] * NO[k,j].conjugate() * NO[k,i];
    return n / (NOoccu.sum());





def derivative(dx, f):
    n = f.size;
    dfdx = np.zeros(n, dtype=np.complex128);
    dfdx[0] = (f[1] - f[n - 2]) / (2 * dx);
    dfdx[n - 1] = dfdx[0];
    for i in range(1, n - 1): dfdx[i] = (f[i + 1] - f[i - 1]) / (2 * dx);
    return dfdx;





def dxFFT(dx, f):
    n = f.size - 1;
    k = 2 * pi * np.fft.fftfreq(n, dx);
    dfdx = np.zeros(f.size, dtype = np.complex128);
    dfdx[:n] = np.fft.fft(f[:n], norm = 'ortho');
    dfdx[:n] = np.fft.ifft(1.0j * k * dfdx[:n], norm = 'ortho');
    dfdx[n] = dfdx[0];
    return dfdx;





def KinE(a2, a1, dx, rho, Orb):
    Morb = rho.shape[0];
    sums = 0.0 + 0.0j;
    for i in range(Morb):

        Ider = dxFFT(dx, Orb[i]);

        for j in range(Morb):
            Jder = dxFFT(dx, Orb[j]);
            sums -= a2 * rho[i,j] * simps(np.conj(Ider) * Jder, dx = dx);
            sums += a1 * rho[i,j] * simps(np.conj(Orb[i]) * Jder, dx = dx);

    return sums;





def PotE(V, dx, rho, Orb):
    Morb = rho.shape[0];
    sums = 0.0 + 0.0j;
    for i in range(Morb):

        for j in range(Morb):
            sums += rho[i,j] * simps(V * np.conj(Orb[i]) * Orb[j], dx = dx);

    return sums;





def IntE(g, dx, rho, Orb):
    Morb = Orb.shape[0];
    M2 = Morb * Morb;
    M3 = Morb * Morb * Morb;
    sums = 0.0 + 0.0j;
    for k in range(Morb):

        for s in range(Morb):

            for q in range(Morb):

                for l in range(Morb):
                    to_int = np.conj(Orb[k] * Orb[s]) * Orb[l] * Orb[q];
                    TwoBody = rho[k + s * Morb + q * M2 + l * M3];
                    sums = sums + TwoBody * simps(to_int, dx = dx);

    return sums * g / 2;





def VonNeumannS(N, RHOeigvals):
    """ (double) = VonNeumannS( # of particles, RHOeigenvalues) """
    return - ( (RHOeigvals.real / N) * np.log(RHOeigvals.real / N) ).sum()
