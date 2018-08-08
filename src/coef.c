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

int main(int argc, char * argv[])
{
    clock_t start, end;
    double cpu_time;

    int N, M;
    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &M);
    
    Carray C = carrDef(NC(N,M));
    Cmatrix rho = cmatDef(M, M);

    int ** NCmat = (int **) malloc((N + 1) * sizeof(int * ));
    for (int i = 0; i < N + 1; i++)
        NCmat[i] = (int *) malloc((M + 1) * sizeof(int));
    for (int i = 0; i < N + 1; i++)
    {
        for (int j = 1; j < M + 1; j++) NCmat[i][j] = NC(i, j);
    }

    carrFill(NC(N,M), (cos(PI / 6) + sin(PI / 6) * I) / sqrt(NC(N, M)), C);

    start = clock();
    for (int i = 0; i < 10; i++) BuildRho(N, M, NCmat, C, rho);
    end = clock();
    cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time to setup rho: %.3fms", cpu_time * 100);


    for (int i = 0; i < N + 1; i++) free(NCmat[i]);
    free(NCmat);
    cmatFree(M, rho);
    free(C);

    printf("\n");
    return 0;
}
