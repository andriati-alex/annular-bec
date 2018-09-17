#include "../include/MCTDHB_configurations.h"





/* ========================================================================
 *
 *       AUXILIAR FUNCTIONS TO CONFIGURE NUMBER OCCUPATION STATES
 *       --------------------------------------------------------
 *
 * A configuration is defined as one of the possibles occupation number 
 * states (Fock vectors).  This is  a  combinatorial problem  on how to
 * fill  "M"  different Single-Particle States (SPS)  with N  available
 * particles.  The routines below  implements a mapping  between  these
 * Fock states and integer numbers, to address coefficients of a  many-
 * body state in Occupation Number Basis (ONB).
 *
 *
 * ======================================================================== */





int fac(int n)
{
    int n_fac = 1;
    for (int i = 1; i < n; i++) n_fac = n_fac * (i + 1);
    return n_fac;
}

int NC(int N, int M)
{
    int n = 1;
    for (int i = N + M - 1; i > N; i --) n = n * i;
    return n / fac(M - 1);
}

int ** MountNCmat(int N, int M)
{
    int i,
        j;

    int ** NCmat = (int ** ) malloc((N + 1) * sizeof(int * ));

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i] = (int * ) malloc((M + 1) * sizeof(int));
    }

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i][0] = 0;
        for (j = 1; j < M + 1; j++) NCmat[i][j] = NC(i, j);
    }

    return NCmat;
}

int ** MountFocks(int N, int M, int ** NCmat)
{
    int k;

    int ** ItoFock = (int **) malloc(NC(N, M) * sizeof(int *));

    for (k = 0; k < NC(N, M); k++)
    {
        ItoFock[k] = (int * ) malloc(M * sizeof(int));
        IndexToFock(k, N, M, NCmat, ItoFock[k]);
    }

    return ItoFock;
}

void IndexToFock(int k, int N, int M, int ** NCmat, int * v)
{
    int x;

    int  i,
         m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    /* -------------------------------------------------------------
     * Put the particles in orbitals while has combinations to spend
    ----------------------------------------------------------------- */
    while ( k > 0 )
    {
        while ( k - NCmat[N][m] < 0 ) m = m - 1;
        x = k - NCmat[N][m];
        while ( x >= 0 )
        {
            v[m] = v[m] + 1; // One more particle in orbital m
            N = N - 1;       // Less one particle to setup
            k = x;
            x = x - NCmat[N][m];
        }
    }

    /* -------------------------------------------------------------
     * Put the particles in orbitals while has combinations to spend
    ----------------------------------------------------------------- */
    for (i = N; i > 0; i--) v[0] = v[0] + 1;
}

int FockToIndex(int N, int M, int ** NCmat, int * v)
{
    int i, n;
    int k = 0;

    /* ---------------------------------------------------
     * Empty one by one orbital starting from the last one
    ------------------------------------------------------ */
    for (i = M - 1; i > 0; i--)
    {
        n = v[i]; // Number of particles in the orbital
        while (n > 0)
        {
            k = k + NCmat[N][i]; // number of combinations needed
            N = N - 1;           // decrease the number of particles
            n = n - 1;
        }
    }

    return k;
}

void JumpMapping(int N, int M, int ** NCmat, int ** IF, int * Map)
{
    int i,
        q,
        k,
        l,
        nc = NCmat[N][M],
        * v = (int *) malloc(M * sizeof(int));

    for (i = 0; i < nc; i++)
    {   // Copy the occupation vector from C[i] coeff.
        for (q = 0; q < M; q++) v[q] = IF[i][q];
        for (k = 0; k < M; k++)
        {   // Take one particle from k state
            if (v[k] < 1) continue;
            for (l = 0; l < M; l++)
            {   // Put one particle in l state
                v[k] -= 1;
                v[l] += 1;
                Map[i + k * nc + l * M * nc] = FockToIndex(N, M, NCmat, v);
                v[k] += 1;
                v[l] -= 1;
            }
        }
    }

    free(v);
}





/* ========================================================================
 *
 *                           <   a*_k   a_l   >
 *                    -------------------------------
 *
 * Once defined a set of Single-Particle Wave Functions (SPWF) a many
 * body state  can be expanded  in a  Occupation Number Configuration
 * Basis (ONCB) whose vector are also named Fock states. The one body
 * density matrix is known as the expected value of 1 creation  and 1
 * annihilation operators for a given many-body state.  Use the basis
 * to express the state and then compute using its coefficients (Cj).
 *
 * ======================================================================== */





void OBrho(int N, int M, int ** NCmat, int ** IF, Carray C, Cmatrix rho)
{
    int i, // int indices to number coeficients
        j,
        k,
        l,
        q,
        nc = NCmat[N][M];

    double mod2;

    double complex RHO;

    int * v; // Occupation vector for each thread

    for (k = 0; k < M; k++)
    {
        RHO = 0;
        #pragma omp parallel shared(k, nc, C, IF) private(i, mod2) \
        reduction(+:RHO)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[i][k];
            }
        }
        rho[k][k] = RHO;

        for (l = k + 1; l < M; l++)
        {
            RHO = 0;
            #pragma omp parallel shared(l,k,M,N,nc,C,NCmat,IF) \
            private(i, j, q, v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));
                #pragma omp for schedule(static)
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 1) continue;
                    for (q = 0; q < M; q++) v[q] = IF[i][q];
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    v[k] += 1;
                    v[l] -= 1;
                    RHO  += conj(C[i]) * C[j] * sqrt((v[l] + 1) * v[k]);
                }
                free(v); // Each thread release its vector
            }
            rho[k][l] = RHO;
        }
    }

    for (k = 0; k < M; k++)
    {   // Use hermiticity to setup lower triangular part
        for (l = k + 1; l < M; l++) rho[l][k] = conj(rho[k][l]);
    }
}





/* ========================================================================
 *
 *                    <   a*_k   a*_s   a_q   a_l   >
 *                    -------------------------------
 *
 * Once defined a set of Single-Particle Wave Functions (SPWF) a many
 * body state  can be expanded  in a  Occupation Number Configuration
 * Basis (ONCB) whose vector are also named Fock states. The two body
 * density matrix is known as the expected value of 2 creation  and 2
 * annihilation operators for a given many-body state.  Use the basis
 * to express the state and then compute using its coefficients (Cj).
 *
 * ======================================================================== */





void TBrho(int N, int M, int ** NCmat, int ** IF, Carray C, Carray rho)
{
    int i, // int indices to number coeficients
        j,
        k,
        s,
        q,
        l,
        t,
        nc = NCmat[N][M];
    
    // Auxiliar to memory access
    int M2 = M * M,
        M3 = M * M * M;

    double mod2,   // |Cj| ^ 2
           sqrtOf; // Factors from the action of creation/annihilation

    double complex RHO;

    int * v; // occupation number vector in each iteration


    /* ---------------------------------------------
     * Rule 1: Creation on k k / Annihilation on k k
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        RHO = 0;
        #pragma omp parallel for private(i,mod2) reduction(+:RHO)
        for (i = 0; i < nc; i++)
        {
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO += mod2 * IF[i][k] * (IF[i][k] - 1);
        }
        rho[k + M * k + M2 * k + M3 * k] = RHO;
    }

    /* ---------------------------------------------
     * Rule 2: Creation on k s / Annihilation on k s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            RHO = 0;
            #pragma omp parallel for private(i,mod2) reduction(+:RHO)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO += mod2 * IF[i][k] * IF[i][s];
            }
            rho[k + s * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + k * M2 + s * M3] = RHO;
            rho[s + k * M + s * M2 + k * M3] = RHO;
            rho[k + s * M + s * M2 + k * M3] = RHO;
        }
    }
    
    /* ---------------------------------------------
     * Rule 3: Creation on k k / Annihilation on q q
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {
            RHO = 0;
            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));
                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = sqrt((v[k] - 1) * v[k] * (v[q] + 1) * (v[q] + 2));
                    v[k] -= 2;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }
                free(v);
            }
            rho[k + k * M + q * M2 + q * M3] = RHO;
        }
    }
 
    /* ---------------------------------------------
     * Rule 4: Creation on k k / Annihilation on k l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            RHO = 0;
            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));
                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 2) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = (v[k] - 1) * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }
                free(v);
            }
            rho[k + k * M + k * M2 + l * M3] = RHO;
        }
    }
    
    /* ---------------------------------------------
     * Rule 5: Creation on k s / Annihilation on s s
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            RHO = 0;
            #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
            {
                v = (int *) malloc(M * sizeof(int));
                #pragma omp for
                for (i = 0; i < nc; i++)
                {
                    if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                    for (t = 0; t < M; t++) v[t] = IF[i][t];
                    sqrtOf = v[s] * sqrt(v[k] * (v[s] + 1));
                    v[k] -= 1;
                    v[s] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    RHO += conj(C[i]) * C[j] * sqrtOf;
                }
                free(v);
            }
            rho[k + s * M + s * M2 + s * M3] = RHO;
        }
    }

    /* -----------------------------------------------------------
     * Rule 6.0: Creation on k k / Annihilation on q l (k < q < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {
                RHO = 0;
                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt( v[k] * (v[k] - 1) * 
                                 (v[q] + 1) * (v[l] + 1) );
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }
                rho[k + k*M + q*M2 + l*M3] = RHO;
            }
        }
    }

    /* -----------------------------------------------------------
     * Rule 6.1: Creation on k k / Annihilation on q l (q < k < l)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (k = q + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {
                RHO = 0;
                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt( v[k] * (v[k] - 1) * 
                                 (v[q] + 1) * (v[l] + 1) );
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }
                rho[k + k*M + q*M2 + l*M3] = RHO;
            }
        }
    }
    
    /* -----------------------------------------------------------
     * Rule 6.2: Creation on k k / Annihilation on q l (q < l < k)
    ------------------------------------------------------------------- */
    for (q = 0; q < M; q++)
    {
        for (l = q + 1; l < M; l++)
        {
            for (k = l + 1; k < M; k++)
            {
                RHO = 0;
                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 2) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = sqrt( v[k] * (v[k] - 1) *
                                 (v[q] + 1) * (v[l] + 1) );
                        v[k] -= 2;
                        v[l] += 1;
                        v[q] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }
                rho[k + k*M + q*M2 + l*M3] = RHO;
            }
        }
    }

    /* -----------------------------------------------------------
     * Rule 7.0: Creation on k s / Annihilation on s l (s < k < l)
    ------------------------------------------------------------------- */
    for (s = 0; s < M; s++)
    {
        for (k = s + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {
                RHO = 0;
                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }
                rho[k + s*M + s*M2 + l*M3] = RHO;
            }
        }
    }

    /* -----------------------------------------------------------
     * Rule 7.1: Creation on k s / Annihilation on s l (k < s < l)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            for (l = s + 1; l < M; l++)
            {
                RHO = 0;
                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }
                rho[k + s*M + s*M2 + l*M3] = RHO;
            }
        }
    }

    /* -----------------------------------------------------------
     * Rule 7.2: Creation on k s / Annihilation on s l (k < l < s)
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            for (s = l + 1; s < M; s++)
            {
                RHO = 0;
                #pragma omp parallel private(i,j,t,sqrtOf,v) reduction(+:RHO)
                {
                    v = (int *) malloc(M * sizeof(int));
                    #pragma omp for
                    for (i = 0; i < nc; i++)
                    {
                        if (IF[i][k] < 1) continue;
                        for (t = 0; t < M; t++) v[t] = IF[i][t];
                        sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                        v[k] -= 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        RHO += conj(C[i]) * C[j] * sqrtOf;
                    }
                    free(v);
                }
                rho[k + s*M + s*M2 + l*M3] = RHO;
            }
        }
    }
        
    /* ---------------------------------------------
     * Rule 8: Creation on k s / Annihilation on q l
    ------------------------------------------------------------------- */
    for (k = 0; k < M; k++)
    {
        for (s = 0; s < M; s++)
        {
            if (s == k) continue;
            for (q = 0; q < M; q++)
            {
                if (q == s || q == k) continue;
                for (l = 0; l < M; l ++)
                {
                    RHO = 0;
                    if (l == k || l == s || l == q) continue;
                    #pragma omp parallel private(i,j,t,sqrtOf,v) \
                    reduction(+:RHO)
                    {
                        v = (int *) malloc(M * sizeof(int));
                        #pragma omp for
                        for (i = 0; i < nc; i++)
                        {
                            if (IF[i][k] < 1 || IF[i][s] < 1) continue;
                            for (t = 0; t < M; t++) v[t] = IF[i][t];
                            sqrtOf = sqrt(v[k]*v[s]*(v[q] + 1)*(v[l] + 1));
                            v[k] -= 1;
                            v[s] -= 1;
                            v[q] += 1;
                            v[l] += 1;
                            j = FockToIndex(N, M, NCmat, v);
                            RHO += conj(C[i]) * C[j] * sqrtOf;
                        }
                        free(v);
                    }
                    rho[k + s*M + q*M2 + l*M3] = RHO;
                } // Finish l loop
            } // Finish q loop
        } // Finish s loop
    } // Finish k loop



    /* ==================================================================== *
     
           TAKE ADVANTAGE OF COMPLEX CONJUGATION AND COMMUTATION RULES
           to fill the rest of elements without call FockToIndex.
      
     * ==================================================================== */



    /****************************** CC Rule 3 ******************************/

    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
            rho[q + q*M + k*M2 + k*M3] = conj(rho[k + k*M + q*M2 + q*M3]);
    }
    
    /****************************** CC Rule 4 ******************************/

    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            rho[k + k*M + l*M2 + k*M3] = rho[k + k*M + k*M2 + l*M3];
            rho[l + k*M + k*M2 + k*M3] = conj(rho[k + k*M + k*M2 + l*M3]);
            rho[k + l*M + k*M2 + k*M3] = rho[l + k*M + k*M2 + k*M3];
        }
    }
    
    /****************************** CC Rule 5 ******************************/

    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            rho[s + k*M + s*M2 + s*M3] = rho[k + s*M + s*M2 + s*M3];
            rho[s + s*M + s*M2 + k*M3] = conj(rho[k + s*M + s*M2 + s*M3]);
            rho[s + s*M + k*M2 + s*M3] = rho[s + s*M + s*M2 + k*M3];
        }
    }
    
    /****************************** CC Rule 6 ******************************/

    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {
                rho[k + k*M + l*M2 + q*M3] = rho[k + k*M + q*M2 + l*M3];
                rho[l + q*M + k*M2 + k*M3] = conj(rho[k + k*M + q*M2 + l*M3]);
                rho[q + l*M + k*M2 + k*M3] = rho[l + q*M + k*M2 + k*M3];
            }
        }
    }

    for (q = 0; q < M; q++)
    {
        for (k = q + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {
                rho[k + k*M + l*M2 + q*M3] = rho[k + k*M + q*M2 + l*M3];
                rho[l + q*M + k*M2 + k*M3] = conj(rho[k + k*M + q*M2 + l*M3]);
                rho[q + l*M + k*M2 + k*M3] = rho[l + q*M + k*M2 + k*M3];
            }
        }
    }
    
    for (q = 0; q < M; q++)
    {
        for (l = q + 1; l < M; l++)
        {
            for (k = l + 1; k < M; k++)
            {
                rho[k + k*M + l*M2 + q*M3] = rho[k + k*M + q*M2 + l*M3];
                rho[l + q*M + k*M2 + k*M3] = conj(rho[k + k*M + q*M2 + l*M3]);
                rho[q + l*M + k*M2 + k*M3] = rho[l + q*M + k*M2 + k*M3];
            }
        }
    }
    
    /****************************** CC Rule 7 ******************************/
        
    for (s = 0; s < M; s++)
    {
        for (k = s + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {
                rho[k + s*M + l*M2 + s*M3] = rho[k + s*M + s*M2 + l*M3];
                rho[s + k*M + l*M2 + s*M3] = rho[k + s*M + s*M2 + l*M3];
                rho[s + k*M + s*M2 + l*M3] = rho[k + s*M + s*M2 + l*M3];
                // First index greater than the last
                rho[l + s*M + s*M2 + k*M3] = conj(rho[k + s*M + s*M2 + l*M3]);
                rho[s + l*M + s*M2 + k*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[s + l*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[l + s*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
            }
        }
    }
    
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            for (l = s + 1; l < M; l++)
            {
                rho[k + s*M + l*M2 + s*M3] = rho[k + s*M + s*M2 + l*M3];
                rho[s + k*M + l*M2 + s*M3] = rho[k + s*M + s*M2 + l*M3];
                rho[s + k*M + s*M2 + l*M3] = rho[k + s*M + s*M2 + l*M3];
                // First index greater than the last
                rho[l + s*M + s*M2 + k*M3] = conj(rho[k + s*M + s*M2 + l*M3]);
                rho[s + l*M + s*M2 + k*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[s + l*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[l + s*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
            }
        }
    }
    
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            for (s = l + 1; s < M; s++)
            {
                rho[k + s*M + l*M2 + s*M3] = rho[k + s*M + s*M2 + l*M3];
                rho[s + k*M + l*M2 + s*M3] = rho[k + s*M + s*M2 + l*M3];
                rho[s + k*M + s*M2 + l*M3] = rho[k + s*M + s*M2 + l*M3];
                // First index greater than the last
                rho[l + s*M + s*M2 + k*M3] = conj(rho[k + s*M + s*M2 + l*M3]);
                rho[s + l*M + s*M2 + k*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[s + l*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[l + s*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
            }
        }
    }

    /* END OF ROUTINE */
}