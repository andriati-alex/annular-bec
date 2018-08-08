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
void IndexToFock(long k, int N, int M, int ** NCmat, int * restrict v);
long FockToIndex(int N, int M, int ** NCmat, int * restrict v);
void BuildRho(int N, int M, int ** NCmat, Carray C, Cmatrix rho);
void BuildRho2(int N, int M, int ** NCmat, Carray C, TwoBodyMat rho);

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

void IndexToFock(long k, int N, int M, int ** NCmat, int * restrict v)
{
    long x;
    int  i,
         m = M - 1;

    for (i = 0; i < M; i++) v[i] = 0;

    /* Put the particles in orbitals while has combinations to spend */
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

    /* Put the rest of particles in the first orbital */
    for (i = N; i > 0; i--) v[0] = v[0] + 1;
}

long FockToIndex(int N, int M, int ** NCmat, int * restrict v)
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

void BuildRho(int N, int M, int ** NCmat, Carray C, Cmatrix rho)
{
    long i, // long indices to number coeficients
         j,
         nc = NC(N, M);
    int  k, // int indices to density matrix (M x M)
         l;

    double mod2,
           sqrtOf;

    /* occupation number vector in each iteration */
    int * restrict v = (int *) malloc(M * sizeof(int));

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

void BuildRho2(int N, int M, int ** NCmat, Carray C, TwoBodyMat rho)
{
    long i, // long indices to number coeficients
         j,
         nc = NC(N, M);
    int  k, // int indices to density matrix (M x M)
         s,
         q,
         l;

    double mod2,    // |Cj| ^ 2
           sqrtOf;  // Factors from the action of creation/annihilation

    /* occupation number vector in each iteration */
    int * restrict v = (int *) malloc(M * sizeof(int));

    for (i = 0; i < nc; i++)
    {
        IndexToFock(i, N, M, NCmat, v);
        mod2 = creal(C[i]) * creal(C[i]) + cimag(C[i]) * cimag(C[i]);

        /****************************** Rule 1 ******************************/
        for (k = 0; k < M; k++) rho[k][k][k][k] += mod2 * v[k] * (v[k] - 1);
        
        /****************************** Rule 2 ******************************/

        for (k = 0; k < M; k++)
        {
            for (s = k + 1; s < M; s++)
            {
                sqrtOf = mod2 * v[k] * v[s];
                rho[k][s][k][s] += sqrtOf;
                rho[s][k][k][s] += sqrtOf;
                rho[s][k][s][k] += sqrtOf;
                rho[k][s][s][k] += sqrtOf;
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
                rho[k][k][q][q] += conj(C[i]) * C[j] * sqrtOf;
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
                rho[k][k][k][l] += conj(C[i]) * C[j] * sqrtOf;
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
                rho[k][s][s][s] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k][k][q][l] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k][k][q][l] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k][k][q][l] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k][s][s][l] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k][s][s][l] += conj(C[i]) * C[j] * sqrtOf;
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
                    rho[k][s][s][l] += conj(C[i]) * C[j] * sqrtOf;
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
                        rho[k][s][q][l] += conj(C[i]) * C[j] * sqrtOf;
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
        for (q = k + 1; q < M; q++) rho[q][q][k][k] = conj(rho[k][k][q][q]);
    }
    
    /****************************** CC Rule 4 ******************************/

    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            rho[k][k][l][k] = rho[k][k][k][l];
            rho[l][k][k][k] = conj(rho[k][k][k][l]);
            rho[k][l][k][k] = rho[l][k][k][k];
        }
    }
    
    /****************************** CC Rule 5 ******************************/

    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            rho[s][k][s][s] = rho[k][s][s][s];
            rho[s][s][s][k] = conj(rho[k][s][s][s]);
            rho[s][s][k][s] = rho[s][s][s][k];
        }
    }
    
    /****************************** CC Rule 6 ******************************/

    for (k = 0; k < M; k++)
    {
        for (q = k + 1; q < M; q++)
        {
            for (l = q + 1; l < M; l++)
            {
                rho[k][k][l][q] = rho[k][k][q][l];
                rho[l][q][k][k] = conj(rho[k][k][q][l]);
                rho[q][l][k][k] = rho[l][q][k][k];
            }
        }
    }
    
    for (q = 0; q < M; q++)
    {
        for (k = q + 1; k < M; k++)
        {
            for (l = k + 1; l < M; l++)
            {
                rho[k][k][l][q] = rho[k][k][q][l];
                rho[l][q][k][k] = conj(rho[k][k][q][l]);
                rho[q][l][k][k] = rho[l][q][k][k];
            }
        }
    }
    
    for (q = 0; q < M; q++)
    {
        for (l = q + 1; l < M; l++)
        {
            for (k = l + 1; k < M; k++)
            {
                rho[k][k][l][q] = rho[k][k][q][l];
                rho[l][q][k][k] = conj(rho[k][k][q][l]);
                rho[q][l][k][k] = rho[l][q][k][k];
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
                rho[k][s][l][s] = rho[k][s][s][l];
                rho[s][k][l][s] = rho[k][s][s][l];
                rho[s][k][s][l] = rho[k][s][s][l];
                rho[l][s][s][k] = conj(rho[k][s][s][l]);
                rho[s][l][s][k] = rho[l][s][s][k];
                rho[s][l][k][s] = rho[l][s][s][k];
                rho[l][s][k][s] = rho[l][s][s][k];
            }
        }
    }
    
    for (k = 0; k < M; k++)
    {
        for (s = k + 1; s < M; s++)
        {
            for (l = s + 1; l < M; l++)
            {
                rho[k][s][l][s] = rho[k][s][s][l];
                rho[s][k][l][s] = rho[k][s][s][l];
                rho[s][k][s][l] = rho[k][s][s][l];
                rho[l][s][s][k] = conj(rho[k][s][s][l]);
                rho[s][l][s][k] = rho[l][s][s][k];
                rho[s][l][k][s] = rho[l][s][s][k];
                rho[l][s][k][s] = rho[l][s][s][k];
            }
        }
    }
    
    for (k = 0; k < M; k++)
    {
        for (l = k + 1; l < M; l++)
        {
            for (s = l + 1; s < M; s++)
            {
                rho[k][s][l][s] = rho[k][s][s][l];
                rho[s][k][l][s] = rho[k][s][s][l];
                rho[s][k][s][l] = rho[k][s][s][l];
                rho[l][s][s][k] = conj(rho[k][s][s][l]);
                rho[s][l][s][k] = rho[l][s][s][k];
                rho[s][l][k][s] = rho[l][s][s][k];
                rho[l][s][k][s] = rho[l][s][s][k];
            }
        }
    }

    free(v);
}

int main(int argc, char * argv[])
{
    clock_t start, end;
    double cpu_time;

    int N, M;
    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &M);
    
    Carray C = carrDef(NC(N,M));
    Cmatrix rho = cmatDef(M, M);
    TwoBodyMat rho2 = TBmatDef(M);

    /* Pre-define values of combinatorial problem to be used */

    int ** NCmat = (int **) malloc((N + 1) * sizeof(int * ));
    for (int i = 0; i < N + 1; i++)
        NCmat[i] = (int *) malloc((M + 1) * sizeof(int));
    for (int i = 0; i < N + 1; i++)
    {
        for (int j = 1; j < M + 1; j++) NCmat[i][j] = NC(i, j);
    }

    /* Array of Coeficients */
    

    carrFill(NC(N,M), (cos(PI / 6) + sin(PI / 6) * I) / sqrt(NC(N, M)), C);

    BuildRho2(N, M, NCmat, C, rho2);

    printf("\n\nrho2[2][1][1][0] = %7.4lf %7.4lf", creal(rho2[2][1][1][0]), cimag(rho2[2][1][1][0]));

    start = clock();
    for (int i = 0; i < 10; i++) BuildRho(N, M, NCmat, C, rho);
    end = clock();
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("\n Time to setup rho: %.3fms", cpu_time * 100);


    for (int i = 0; i < N + 1; i++) free(NCmat[i]);
    free(NCmat);
    cmatFree(M, rho);
    free(C);

    printf("\n");
    return 0;
}
