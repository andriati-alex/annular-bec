#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "../include/array_memory.h"
#include "../include/array_operations.h"

#define PI 3.141592653589793

long fac(int n);
long NC(int N, int M);
void IndexToFock(long k, int N, int M, long ** NCmat, int * v);
long FockToIndex(int N, int M, long ** NCmat, int * v);
void BuildRho(int N, int M, long ** NCmat, Carray C, Cmatrix rho);
void BuildRho2(int N, int M, long ** NCmat, Carray C, Carray rho);

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

void IndexToFock(long k, int N, int M, long ** NCmat, int * v)
{
    long x;
    int  i,
         m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    /*** Put the particles in orbitals while has combinations to spend ***/

    while ( k > 0 )
    {
        while ( k - NCmat[N][m] < 0 ) m = m - 1;
        x = k - NCmat[N][m];
        while ( x >= 0 )
        {
            v[m] = v[m] + 1;
            N = N - 1;
            k = x;
            x = x - NCmat[N][m];
        }
    }

    /***   Put the rest of particles in the first orbital   ***/

    for (i = N; i > 0; i--) v[0] = v[0] + 1;
}

long FockToIndex(int N, int M, long ** NCmat, int * v)
{
    int  i, n;
    long k = 0;

    /* Empty one by one orbital starting from the last one */
    for (i = M - 1; i > 0; i--)
    {
        n = v[i]; // Number of particles in current orbital
        while (n > 0)
        {
            k = k + NCmat[N][i]; // number of combinations needed
            N = N - 1;           // decrease the number of particles
            n = n - 1;
        }
    }
    return k;
}

void BuildRho(int N, int M, long ** NCmat, Carray C, Cmatrix rho)
{
    long i, // long indices to number coeficients
         j,
         nc = NCmat[N][M];
    int  k, // int indices to density matrix (M x M)
         l;

    double mod2,
           sqrtOf;

    /* occupation number vector in each iteration */
    int * v = (int *) malloc(M * sizeof(int));

    /* Initialize the density matrix with zero */
    for (k = 0; k < M; k++)
    {
        rho[k][k] = 0;
        for (l = k + 1; l < M; l++) rho[k][l] = 0;
    }

    for (i = 0; i < nc; i++)
    {
        IndexToFock(i, N, M, NCmat, v);
        mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
        for (k = 0; k < M; k++)
        {
            rho[k][k] += v[k] * mod2;
            for (l = k + 1; l < M; l++)
            {
                if (v[k] < 1) continue; // annihilation on vaccumm
                sqrtOf = sqrt( (v[l] + 1) * v[k] );
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);
                rho[k][l] += conj(C[i]) * C[j] * sqrtOf;
                v[k] += 1;
                v[l] -= 1;
            }
        }
    }

    /* Use hermiticity of density matrix to fill lower part */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++) rho[l][k] = conj(rho[k][l]);
    }

    free(v);
}

void BuildRhoP(int N, int M, long ** NCmat, Carray C, Cmatrix rho)
{
    long i, // long indices to number coeficients
         j,
         nc = NCmat[N][M];
    int  k, // int indices to density matrix (M x M)
         l;

    double mod2,
           sqrtOf;

    double complex RHO_kk,
                   RHO_kl;

    /* occupation number vector in each iteration */
    int * v = (int *) malloc(M * sizeof(int));

    for (k = 0; k < M; k++)
    {
        RHO_kk = 0;
        for (i = 0; i < nc; i++)
        {
            IndexToFock(i, N, M, NCmat, v);
            mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);
            RHO_kk += mod2 * v[k];
        }
        rho[k][k] = RHO_kk;

        for (l = k + 1; l < M; l++)
        {
            RHO_kl = 0;
            for (i = 0; i < nc; i++)
            {
                IndexToFock(i, N, M, NCmat, v);
                if (v[k] < 1) continue;
                sqrtOf = sqrt((v[l] + 1) * v[k]);
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);
                RHO_kl += conj(C[i]) * C[j] * sqrtOf;
                v[k] += 1;
                v[l] -= 1;
            }
            rho[k][l] = RHO_kl;
        }
    }

    /* Use hermiticity of density matrix to fill lower part */
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++) rho[l][k] = conj(rho[k][l]);
    }

    free(v);
}

void BuildRho2(int N, int M, long ** NCmat, Carray C, Carray rho)
{
    long i, // long indices to number coeficients
         j,
         nc = NCmat[N][M];
    int  k, // int indices to density matrix (M x M)
         s,
         q,
         l;

    int M2 = M * M,
        M3 = M * M * M;

    double mod2,    // |Cj| ^ 2
           sqrtOf;  // Factors from the action of creation/annihilation

    /***   occupation number vector in each iteration   ***/

    int * v = (int *) malloc(M * sizeof(int));

    // Clean up any value before start
    for (k = 0; k < M * M * M * M; k++) rho[k] = 0;

    for (i = 0; i < nc; i++)
    {
        IndexToFock(i, N, M, NCmat, v);
        mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);

        /****************************** Rule 1 ******************************/
        for (k = 0; k < M; k++)
            rho[k + M * k + M2 * k + M3 * k] += mod2 * v[k] * (v[k] - 1);
        
        /****************************** Rule 2 ******************************/

        for (k = 0; k < M; k++)
        {
            for (s = k + 1; s < M; s++)
            {
                sqrtOf = mod2 * v[k] * v[s];
                rho[k + s * M + k * M2 + s * M3] += sqrtOf;
                rho[s + k * M + k * M2 + s * M3] += sqrtOf;
                rho[s + k * M + s * M2 + k * M3] += sqrtOf;
                rho[k + s * M + s * M2 + k * M3] += sqrtOf;
            }
        }

        /****************************** Rule 3 ******************************/

        for (k = 0; k < M; k++)
        {
            if (v[k] < 2) continue;
            for (q = k + 1; q < M; q++)
            {
                sqrtOf = sqrt((v[k] - 1) * v[k] * (v[q] + 1) * (v[q] + 2));
                v[k] -= 2;
                v[q] += 2;
                j = FockToIndex(N, M, NCmat, v);
                rho[k + k * M + q * M2 + q * M3] += conj(C[i]) * C[j] * sqrtOf;
                v[k] += 2;
                v[q] -= 2;
            }
        }
        
        /****************************** Rule 4 ******************************/

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (l = k + 1; l < M; l++)
            {
                sqrtOf = (v[k] - 1) * sqrt(v[k] * (v[l] + 1));
                v[k] -= 1;
                v[l] += 1;
                j = FockToIndex(N, M, NCmat, v);
                rho[k + k * M + k * M2 + l * M3] += conj(C[i]) * C[j] * sqrtOf;
                v[k] += 1;
                v[l] -= 1;
            }
        }
        
        /****************************** Rule 5 ******************************/
        
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                sqrtOf = v[s] * sqrt(v[k] * (v[s] + 1));
                v[k] -= 1;
                v[s] += 1;
                j = FockToIndex(N, M, NCmat, v);
                rho[k + s * M + s * M2 + s * M3] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k + k*M + q*M2 + l*M3] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k + k*M + q*M2 + l*M3] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k + k*M + q*M2 + l*M3] += conj(C[i]) * C[j] * sqrtOf;
                    v[k] += 2;
                    v[l] -= 1;
                    v[q] -= 1;
                }
            }
        }

        /****************************** Rule 7 ******************************/

        for (s = 0; s < M; s++)
        {
            for (k = s + 1; k < M; k++)
            {
                if (v[k] < 1) continue;
                for (l = k + 1; l < M; l++)
                {
                    sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    rho[k + s*M + s*M2 + l*M3] += conj(C[i]) * C[j] * sqrtOf;
                    v[k] += 1;
                    v[l] -= 1;
                }
            }
        }

        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (s = k + 1; s < M; s++)
            {
                for (l = s + 1; l < M; l++)
                {
                    sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    rho[k + s*M + s*M2 + l*M3] += conj(C[i]) * C[j] * sqrtOf;
                    v[k] += 1;
                    v[l] -= 1;
                }
            }
        }
        
        for (k = 0; k < M; k++)
        {
            if (v[k] < 1) continue;
            for (l = k + 1; l < M; l++)
            {
                for (s = l + 1; s < M; s++)
                {
                    sqrtOf = v[s] * sqrt(v[k] * (v[l] + 1));
                    v[k] -= 1;
                    v[l] += 1;
                    j = FockToIndex(N, M, NCmat, v);
                    rho[k + s*M + s*M2 + l*M3] += conj(C[i]) * C[j] * sqrtOf;
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
                        sqrtOf = sqrt(v[k]*v[s]*(v[q] + 1)*(v[l] + 1));
                        v[k] -= 1;
                        v[s] -= 1;
                        v[q] += 1;
                        v[l] += 1;
                        j = FockToIndex(N, M, NCmat, v);
                        rho[k+s*M+q*M2+l*M3] += conj(C[i]) * C[j] * sqrtOf;
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

        /* Go to next coeficient in the sum */
    }

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

    free(v);
}

int main(int argc, char * argv[])
{
    clock_t start, end;
    double cpu_time;

    int Npar, // Number of particles
        Morb, // Number of orbitals
        i,
        j;

    sscanf(argv[1], "%d", &Npar);
    sscanf(argv[2], "%d", &Morb);

    /***   Pre-define values of combinatorial problem to be used   ***/

    printf("\nTotal number of coeficients: %ld\n", NC(Npar, Morb));

    long ** NCmat;
    NCmat = (long ** ) malloc((Npar + 1) * sizeof(long * ));

    for (i = 0; i < Npar + 1; i++)
        NCmat[i] = (long * ) malloc((Morb + 1) * sizeof(long));

    for (i = 0; i < Npar + 1; i++)
    {
        NCmat[i][0] = 0;
        for (j = 1; j < Morb + 1; j++) NCmat[i][j] = NC(i, j);
    }

    Carray C = carrDef( NC( Npar, Morb ) );     // Coeficients of superposition
    Cmatrix rho = cmatDef( Morb,  Morb);        // onde-body density matrix
    Carray rho2 = carrDef(Morb*Morb*Morb*Morb); // two-body density matrix

    /***   Setup array of coeficients   ***/

    carrFill(NC(Npar, Morb), 
             (cos(PI / 6) + sin(PI / 6) * I) / sqrt(NC(Npar, Morb)), C);

    /***   Time performance to setup One-Body part   ***/

    start = clock();

    for (i = 0; i < 10; i++) BuildRhoP(Npar, Morb, NCmat, C, rho);

    end = clock();

    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC * 100;

    printf("\n\tTime to setup rho  (Npar = %d, Morb = %d): %.3fms", 
            Npar, Morb, cpu_time);

    /***   Time performance to setup Two-Body part   ***/

    start = clock();

    BuildRho2(Npar, Morb, NCmat, C, rho2);
    BuildRho2(Npar, Morb, NCmat, C, rho2);

    end = clock();
    
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC * 500;

    printf("\n\n\tTime to setup rho2 (Npar = %d, Morb = %d): %.3fms", 
            Npar, Morb, cpu_time);

    /***   Release Memory   ***/

    for (i = 0; i < Npar + 1; i++) free(NCmat[i]); free(NCmat);

    cmatFree(Morb, rho);

    free(rho2);

    free(C);

    /***   Print Some Value for the case N = M = 3   ***/

    Npar = 3;
    Morb = 3;
    
    NCmat = (long ** ) malloc((Npar + 1) * sizeof(long * ));

    for (i = 0; i < Npar + 1; i++)
        NCmat[i] = (long * ) malloc((Morb + 1) * sizeof(long));

    for (i = 0; i < Npar + 1; i++)
    {
        NCmat[i][0] = 0;
        for (j = 1; j < Morb + 1; j++) NCmat[i][j] = NC(i, j);
    }
    
    C = carrDef( NC( Npar, Morb ) );

    carrFill(NC(Npar, Morb), 1, C);

    C[2] = - 1 * I; C[5] = 1 + 1 * I; C[6] = 1 - 1 * I; C[8] = 0;
    
    rho = cmatDef( Morb,  Morb);         // onde-body density matrix
    rho2 = carrDef(Morb*Morb*Morb*Morb); // two-body density matrix
    
    BuildRhoP(Npar, Morb, NCmat, C, rho);
    BuildRhoP(Npar, Morb, NCmat, C, rho);
    
    BuildRho2(Npar, Morb, NCmat, C, rho2);
    BuildRho2(Npar, Morb, NCmat, C, rho2);

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

    
    /***   Release Memory   ***/

    for (i = 0; i < Npar + 1; i++) free(NCmat[i]); free(NCmat);

    cmatFree(Morb, rho);

    free(rho2);

    free(C);

    printf("\n\n");
    return 0;
}
