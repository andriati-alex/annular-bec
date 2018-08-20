#include "../include/coef_routines.h"

long fac(int n)
{
    long n_fac = 1;
    for (int i = 1; i < n; i++) n_fac = n_fac * (i + 1);
    return n_fac;
}

long NC(int N, int M)
{
    long n = 1;
    for (int i = N + M - 1; i > N; i --) n = n * i;
    return n / fac(M - 1);
}

long ** MountNCmat(int N, int M)
{
    int i,
        j;

    long ** NCmat = (long ** ) malloc((N + 1) * sizeof(long * ));

    for (i = 0; i < N + 1; i++)
        NCmat[i] = (long * ) malloc((M + 1) * sizeof(long));

    for (i = 0; i < N + 1; i++)
    {
        NCmat[i][0] = 0;
        for (j = 1; j < M + 1; j++) NCmat[i][j] = NC(i, j);
    }

    return NCmat;
}

int ** MountFocks(int N, int M, long ** NCmat)
{
    long k;

    int ** ItoFock = (int **) malloc(NC(N, M) * sizeof(int *));

    for (k = 0; k < NC(N, M); k++)
    {
        ItoFock[k] = (int * ) malloc(M * sizeof(int));
        IndexToFock(k, N, M, NCmat, ItoFock[k]);
    }

    return ItoFock;
}

void IndexToFock(long k, int N, int M, long ** NCmat, int * v)
{
    long x;

    int  i,
         m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    /***  Put the particles in orbitals while has combinations to spend  ***/

    while ( k > 0 )
    {
        while ( k - NCmat[N][m] < 0 ) m = m - 1;
        x = k - NCmat[N][m];
        while ( x >= 0 )
        {
            v[m] = v[m] + 1; // One more particle in orbital m
            N = N - 1;       // Less one particle still to setup
            k = x;
            x = x - NCmat[N][m];
        }
    }

    /***         Put the rest of particles in the first orbital         ***/

    for (i = N; i > 0; i--) v[0] = v[0] + 1;
}

long FockToIndex(int N, int M, long ** NCmat, int * v)
{
    int  i, n;
    long k = 0;

    /*******  Empty one by one orbital starting from the last one  *******/

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

void JumpMapping(int N, int M, long ** NCmat, int ** IF, long * Map)
{
    long i,
         nc = NCmat[N][M];

    int  q,
         k,
         l;

    int * v = (int *) malloc(M * sizeof(int));

    for (i = 0; i < nc; i++)
    {
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (l = 0; l < M; l++)
            {
                for (q = 0; q < M; q++) v[q] = IF[i][q];
                v[k] -= 1;
                v[l] += 1;
                Map[i + k * nc + l * M * nc] = FockToIndex(N, M, NCmat, v);
            }
        }
    }

    free(v);
}


void OBrho(int N, int M, long ** NCmat, int ** IF, Carray C, Cmatrix rho)
{
    long i, // long indices to number coeficients
         j,
         nc = NCmat[N][M];
    int  k, // int indices to density matrix (M x M)
         l,
         q;

    double mod2,
           sqrtOf;

    double complex RHO_kk,
                   RHO_kl;

    int * v;

    for (k = 0; k < M; k++)
    {
        RHO_kk = 0;
        #pragma omp parallel shared(k,nc,C,IF) private(i,mod2) \
        reduction(+:RHO_kk)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < nc; i++)
            {
                mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
                RHO_kk += mod2 * IF[i][k];
            }
        }
        rho[k][k] = RHO_kk;

        for (l = k + 1; l < M; l++)
        {
            RHO_kl = 0;
            #pragma omp parallel shared(l,k,M,N,nc,C,NCmat,IF) \
            private(i,j,q,v) reduction(+:RHO_kl)
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
                    RHO_kl += conj(C[i]) * C[j] * sqrt((v[l] + 1) * v[k]);
                }
                free(v); // Each thread release its vector
            }
            rho[k][l] = RHO_kl;
        }
    }

    /*******  Use hermiticity of density matrix to fill lower part  *******/

    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++) rho[l][k] = conj(rho[k][l]);
    }

}

void TBrho(int N, int M, long ** NCmat, int ** IF, Carray C, Carray rho)
{
    long i, // long indices to number coeficients
         j,
         nc = NCmat[N][M];
    int  k, // int indices to density matrix (M x M x M x M)
         s,
         q,
         l,
         t; // To copy the fock vector

    int M2 = M * M,
        M3 = M * M * M;

    double mod2,   // |Cj| ^ 2
           sqrtOf; // Factors from the action of creation/annihilation

    double complex RHO;

    int * v; // occupation number vector in each iteration


    /****************************** Rule 1 ******************************/

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

    /****************************** Rule 2 ******************************/

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
    
    /****************************** Rule 3 ******************************/

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
        
    /****************************** Rule 4 ******************************/

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
                    if (IF[i][k] < 1) continue;
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
    
    /****************************** Rule 5 ******************************/

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
                    if (IF[i][k] < 1) continue;
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

    /****************************** Rule 6 ******************************/
        
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

    /****************************** Rule 7 ******************************/

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
        
    /****************************** Rule 8 ******************************/
 
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
                    rho[k+s*M+q*M2+l*M3] = RHO;
                } // Finish l loop
            } // Finish q loop
        } // Finish s loop
    } // Finish k loop

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
                rho[l + s*M + s*M2 + k*M3] = conj(rho[k + s*M + s*M2 + l*M3]);
                rho[s + l*M + s*M2 + k*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[s + l*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
                rho[l + s*M + k*M2 + s*M3] = rho[l + s*M + s*M2 + k*M3];
            }
        }
    }

    // END ROUTINE
}

void RHSofODES(int N, int M, int ** IF, long ** NCmat,
               Cmatrix Ho, Carray Hint, Carray C, Carray rhs)
{
    long i, // index of coeficient
         j;

    int  k, // enumerate orbitals
         l,
         s,
         q,
         t;

    int  M2 = M * M,
         M3 = M * M * M;

    int * v = (int * ) malloc(M * sizeof(int));

    double sqrtOf;

    double complex rhsI; // line step value of right-hand-side

    for (i = 0; i < NCmat[N][M]; i++)
    {
        rhsI = 0;

        for (t = 0; t < M; t++) v[t] = IF[i][t];

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            rhsI += Ho[k][k] * v[k] * C[i];
            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = sqrt(v[k] * (v[l] + 1));
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);
                rhsI += Ho[k][l] * sqrtOf * C[j];
                v[k] += 1;
                v[l] -= 1;
            }
        }

        /****************************** Rule 1 ******************************/
        for (k = 0; k < M; k++)
            rhsI += Hint[k + M*k + M2*k + M3*k] * C[i] * v[k] * (v[k] - 1);
        
        /****************************** Rule 2 ******************************/

        for (k = 0; k < M; k++)
        {
            for (s = k + 1; s < M; s++)
            {
                sqrtOf = v[k] * v[s];
                rhsI += Hint[k + s*M + k*M2 + s*M3] * sqrtOf * C[i];
                rhsI += Hint[s + k*M + k*M2 + s*M3] * sqrtOf * C[i];
                rhsI += Hint[s + k*M + s*M2 + k*M3] * sqrtOf * C[i];
                rhsI += Hint[k + s*M + s*M2 + k*M3] * sqrtOf * C[i];
            }
        }

        /****************************** Rule 3 ******************************/

        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = 0; q < M; q++)
            {
                if (q == k) continue;
                sqrtOf = sqrt((v[k] - 1) * v[k] * (v[q] + 1) * (v[q] + 2));
                v[k] -= 2;
                v[q] += 2;
                j = FockToIndex(N, M, NCmat, v);
                rhsI += Hint[k + k * M + q * M2 + q * M3] * C[j] * sqrtOf;
                v[k] += 2;
                v[q] -= 2;
            }
        }

        /****************************** Rule 4 ******************************/

        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (l = 0; l < M; l++)
            {
                if (l == k) continue;
                sqrtOf = (v[k] - 1) * sqrt(v[k] * (v[l] + 1));
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);
                rhsI += Hint[k + k * M + k * M2 + l * M3] * C[j] * sqrtOf;
                rhsI += Hint[k + k * M + l * M2 + k * M3] * C[j] * sqrtOf;
                v[k] += 1;
                v[l] -= 1;
            }
        }

        /****************************** Rule 5 ******************************/

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = 0; s < M; s++)
            {
                if (s == k) continue;
                sqrtOf = v[s] * sqrt(v[k] * (v[s] + 1));
                v[k] -= 1;
                v[s] += 1;
                j = FockToIndex(N, M, NCmat, v);
                rhsI += Hint[k + s * M + s * M2 + s * M3] * C[j] * sqrtOf;
                rhsI += Hint[s + k * M + s * M2 + s * M3] * C[j] * sqrtOf;
                v[k] += 1;
                v[s] -= 1;
            }
        }

        /****************************** Rule 6 ******************************/

        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = k + 1; q < M; q++)
            {
                for (l = q + 1; l < M; l++)
                {
                    sqrtOf = sqrt(v[k] * (v[k] - 1) * (v[q] + 1) * (v[l] + 1));
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    rhsI += Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    rhsI += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }
        }

        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 2) continue;
                for (l = k + 1; l < M; l++)
                {
                    sqrtOf = sqrt(v[k] * (v[k] - 1) * (v[q] + 1) * (v[l] + 1));
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    rhsI += Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    rhsI += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }
        }

        for (q = 0; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {
                for (k = l + 1; k < M; k++)
                {
                    if (v[k] < 2) continue;
                    sqrtOf = sqrt(v[k] * (v[k] - 1) * (v[q] + 1) * (v[l] + 1));
                    v[k] -= 2;
                    v[l] += 1;
                    v[q] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    rhsI += Hint[k + k*M + q*M2 + l*M3] * C[j] * sqrtOf;
                    rhsI += Hint[k + k*M + l*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }
        }
        
        /****************************** Rule 6.1 ******************************/

        for (q = 0; q < M; q++)
        {
            for (k = q + 1; k < M; k++)
            {
                if (v[k] < 1) continue;
                for (s = k + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt(v[k] * v[s] * (v[q] + 1) * (v[q] + 2));
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    rhsI += Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    rhsI += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (q = k + 1; q < M; q++)
            {
                for (s = q + 1; s < M; s++)
                {
                    if (v[s] < 1) continue;
                    sqrtOf = sqrt(v[k] * v[s] * (v[q] + 1) * (v[q] + 2));
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    rhsI += Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    rhsI += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                if (v[s] < 1) continue;
                for (q = s + 1; q < M; q++)
                {
                    sqrtOf = sqrt(v[k] * v[s] * (v[q] + 1) * (v[q] + 2));
                    v[k] -= 1;
                    v[s] -= 1;
                    v[q] += 2;
                    j = FockToIndex(N, M, NCmat, v);
                    rhsI += Hint[k + s*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    rhsI += Hint[s + k*M + q*M2 + q*M3] * C[j] * sqrtOf;
                    v[k] += 1;
                    v[s] += 1;
                    v[q] -= 2;
                }
            }
        }

        /****************************** Rule 7 ******************************/

        for (s = 0; s < M; s++)
        {
            // if (v[s] < 1) continue // may improve performance 
            for (k = 0; k < M; k++)
            {
                if (v[k] < 1 || k == s) continue;
                for (l = k + 1; l < M; l++)
                {
                    if (l == k || l == s) continue;
                    sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    rhsI += Hint[k + s*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    rhsI += Hint[s + k*M + s*M2 + l*M3] * C[j] * sqrtOf;
                    rhsI += Hint[s + k*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    rhsI += Hint[k + s*M + l*M2 + s*M3] * C[j] * sqrtOf;
                    v[k] += 1;
                    v[l] -= 1;
                }
            }
        }
        
        /****************************** Rule 8 ******************************/

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = 0; s < M; s++)
            {
                if (v[s] < 1 || s == k) continue;
                for (q = 0; q < M; q++)
                {
                    if (q == s || q == k) continue;
                    for (l = 0; l < M; l ++)
                    {
                        if (l == k || l == s || l == q) continue;
                        sqrtOf = sqrt(v[k] * v[s] * (v[q] + 1) * (v[l] + 1));
                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        rhsI += Hint[k + s*M + q*M2 + l*M3] * C[j] * sqrtOf;
                        v[k] += 1;
                        v[s] += 1;
                        v[q] -= 1;
                        v[l] -= 1;
                    }
                    // Finish l
                }
                // Finish q
            }
            // Finish s
        }
        // Finish k

        rhs[i] = rhsI;
    }

    free(v);
}

int main(int argc, char * argv[])
{
    omp_set_num_threads(omp_get_max_threads() / 2);

    clock_t start, end;
    double time_used;
    double startP;

    int Npar, // Number of particles
        Morb, // Number of orbitals
        i,
        j;

    long k;

    sscanf(argv[1], "%d", &Npar);
    sscanf(argv[2], "%d", &Morb);



    /***   Pre-define values of combinatorial problem to be used   ***/



    printf("\nTotal number of coeficients: %ld\n", NC(Npar, Morb));

    long ** NCmat = MountNCmat(Npar, Morb);



    /*********                 Main Variables                 *********/



    Carray C = carrDef( NC( Npar, Morb ) );     // Coeficients of superposition
    Cmatrix rho = cmatDef( Morb,  Morb);        // onde-body density matrix
    Carray rho2 = carrDef(Morb*Morb*Morb*Morb); // two-body density matrix

    /***   Setup array of coeficients   ***/

    for (k = 0; k < NCmat[Npar][Morb]; k++)
    {
        C[k] = cos(6 * PI / NCmat[Npar][Morb]) + \
               sin(6 * PI / NCmat[Npar][Morb]) * I;
        C[k] = C[k] / NCmat[Npar][Morb]; // Normalize
    }



    /***         All ocuppation number vectors in a Matrix          ***/



    int ** ItoFock = MountFocks(Npar, Morb, NCmat);



    /***          Time performance to setup One-Body part          ***/



    startP = omp_get_wtime();

    for (i = 0; i < 10; i++) OBrho(Npar, Morb, NCmat, ItoFock, C, rho);

    time_used = (double) (omp_get_wtime() - startP) * 100;

    printf("\n\tTime to setup rho  (Npar = %d, Morb = %d): %.3fms", 
            Npar, Morb, time_used);



    /***         Time performance to setup Two-Body part         ***/



    startP = omp_get_wtime();

    TBrho(Npar, Morb, NCmat, ItoFock, C, rho2);
    TBrho(Npar, Morb, NCmat, ItoFock, C, rho2);

    time_used = (double) (omp_get_wtime() - startP) * 500;

    printf("\n\n\tTime to setup rho2 (Npar = %d, Morb = %d): %.3fms", 
            Npar, Morb, time_used);



    /***                    Release Memory                    ***/



    for (i = 0; i < Npar + 1; i++) free(NCmat[i]); free(NCmat);
    
    for (k = 0; k < NC(Npar, Morb); k++) free(ItoFock[k]); free(ItoFock);

    cmatFree(Morb, rho);

    free(rho2);

    free(C);



    /***       Print Some Values for the case N = M = 3       ***/



    Npar = 3;
    Morb = 3;
    
    NCmat = MountNCmat(Npar, Morb);
    
    ItoFock = MountFocks(Npar, Morb, NCmat);
    
    C = carrDef( NC( Npar, Morb ) );

    carrFill(NC(Npar, Morb), 1, C);

    C[2] = - 1 * I; C[5] = 1 + 1 * I; C[6] = 1 - 1 * I; C[8] = 0;
    
    rho = cmatDef( Morb, Morb );                // onde-body density matrix
    rho2 = carrDef(Morb * Morb * Morb * Morb);  // two-body density matrix
    
    OBrho(Npar, Morb, NCmat, ItoFock, C, rho);
    OBrho(Npar, Morb, NCmat, ItoFock, C, rho);

    TBrho(Npar, Morb, NCmat, ItoFock, C, rho2);
    TBrho(Npar, Morb, NCmat, ItoFock, C, rho2);

    printf("\n\nMatrix rho for N = M = 3\n");
    for (i = 0; i < Morb; i++)
    {
        printf("\n\t|");
        for (j = 0; j < Morb; j++)
        {
            printf("  %7.3lf %7.3lfj |", creal(rho[i][j]), cimag(rho[i][j]));
        }
    }

    printf("\n\nrho2[1,1,0,2] = %7.3lf %7.3lf", creal(rho2[1 + 3 + 2*27]),
            cimag(rho2[1 + 3 + 2*27]));
    
    printf("\n\nrho2[1,2,2,0] = %7.3lf %7.3lf", creal(rho2[1 + 6 + 18]),
            cimag(rho2[1 + 6 + 18]));
    
    printf("\n\nrho2[1,2,2,1] = %7.3lf %7.3lf", creal(rho2[1 + 6 + 18 + 27]),
            cimag(rho2[1 + 6 + 18 + 27]));
    
    printf("\n\nrho2[1,1,1,0] = %7.3lf %7.3lf", creal(rho2[1 + 3 + 9]),
            cimag(rho2[1 + 3 + 9]));



    /**************************   Release Memory   **************************/



    for (i = 0; i < Npar + 1; i++) free(NCmat[i]); free(NCmat);
    
    for (k = 0; k < NC(Npar, Morb); k++) free(ItoFock[k]); free(ItoFock);

    cmatFree(Morb, rho);

    free(rho2);

    free(C);

    printf("\n\n");
    return 0;
}
