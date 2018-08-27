#include "../include/MCTDHB_integrator.h"

struct orbitals
{
    int Morb;   // # of orbitals
    int Mpos;   // # discretized positions (# of divisions + 1)
    double dx;  // spatial step
    double xi;
    double xf;
    Cmatrix O;  // Values of orbitals on each position (Morb x M)
};

typedef struct orbitals * MCorital;

struct coeficients
{
    int Npar;       // # of particles
    int Morb;       // # of orbitals
    long nc;        // all possible fock states
    int ** IF;      // All occu. vector by row
    long ** NCmat;  // All outcomes from combinatorial
    Carray C;       // Array of coeficients
};

typedef struct coeficients * MCcoef;

struct eq_parameters
{
    double a2;          // Second order derivative term
    double inter;       // interaction 'strength'
    double complex a1;  // First order derivative term
    Rarray V;           // Values of potential in each position
};

typedef struct eq_parameters * EqSetup;

MCorbital AllocOrb(int M, int Mpos, double xi, double xf)
{
    MCorbital psi = (MCorbital) malloc(sizeof(struct orbitals));
    psi->Morb = M;
    psi->Mpos = Mpos;
    psi->xi = xi;
    psi->xf = xf;
    psi->dx = (xf - xi) / (Mpos - 1);
    psi->O = cmatDef(M, Mpos);
    return psi;
}

MCcoef AllocCoef(int Npar, int Morb)
{
    MCcoef Coef = (MCcoef) malloc(sizeof(struct coeficients));
    Coef->Npar = Npar;
    Coef->Morb = Morb;
    Coef->NCmat = MountNCmat(Npar, Morb);
    Coef->nc = Coef->NCmat[Npar][Morb];
    Coef->IF = MountFocks(Npar, Morb, Coef->NCmat);
    Coef->C = carrDef(Coef->nc);
    return Coef;
}

EqSetup AllocEq(double a2, double complex a1, double inter, Rarray V)
{
    EqSetup eq = (EqSetup) malloc(sizeof(struct eq_parameters));
    eq->a2 = a2;
    eq->a1 = a1;
    eq->inter = inter;
    eq->V = V;
    return eq;
}

void SetupHo(int Morb, int Mpos, Cmatrix Omat, double dx,
             double a2, double complex a1, Rarray V, Cmatrix Ho)
{
    int i,
        j,
        k;

    Carray ddxi  = carrDef(Mpos);
    Carray ddxj  = carrDef(Mpos);
    Carray toInt = carrDef(Mpos);

    for (i = 0; i < Morb; i++)
    {
        dxCyclic(Mpos, Omat[i], dx, ddxi);
        for (j = i; j < Morb; j++)
        {
            dxCyclic(Mpos, Omat[j], dx, ddxj);
            for (k = 0; k < Mpos; k++)
            {
                toInt[k] = - a2 * conj(ddxi[k]) * ddxj[k]    \
                           + a1 * conj(Omat[i][k]) * ddxj[k] \
                           + V[k] * conj(Omat[i][k]) * Omat[j][k];
            }
            Ho[i][j] = Csimps(Mpos, toInt, dx);
            Ho[j][i] = conj(H[i][j]);
        }
    }

    free(ddxi); free(ddxj); free(toInt);
}

void SetupHint(int Morb, int Mpos, Cmatrix Omat, double dx,
               double inter, Cmatrix Hint)
{
    int i,
        k,
        s,
        q,
        l;

    double complex Integral;

    Carray toInt = carrDef(N);

    for (k = 0; k < Morb; k++)
    {
        for (s = k; s < Morb; s++)
        {
            for (q = 0; q < Morb; q++)
            {
                for (l = q; l < Morb; l++)
                {
                    for (i = 0; i < Mpos; i++)
                    {
                        toInt[i] = conj(Omat[k][i] * Omat[s][i]) * \
                                   Omat[q][i] * Omag[l][i];
                    }
                    Integral = inter * Csimps(Mpos, toInt, dx);
                    Hint[k + s * M + q * M2 + l * M3] = Integral;
                    Hint[k + s * M + l * M2 + q * M3] = Integral;
                    Hint[s + k * M + q * M2 + l * M3] = Integral;
                    Hint[s + k * M + l * M2 + q * M3] = Integral;
                }   // Take advantage of the symmetry k <--> s
            }
        }
    }

    free(toInt);
}

double complex Proj_Hint(int M, int k, int i,
                         Cmatrix rho_inv, Carray rho2, Carray Hint)
{
    int j,
        s,
        q,
        l,
        i_rho2,
        i_Hint;

    double complex ans = 0;

    for (j = 0; j <M; j++)
    {
        for (s = 0; s < M; s++)
        {
            for (q = 0; q < M; q++)
            {
                for (l = 0; l < M; l++)
                {
                    i_rho2 = j + s * M + q * M * M + l * M * M * M;
                    i_Hint = i + s * M + l * M * M + q * M * M * M;
                    ans += rho_inv[k][j] * rho2[i_rho2] * Hint[i_Hint];
                }
            }
        }
    }

    return ans;
}

double complex NonLinear(int M, int k, int n,
                         Cmatrix Omat, Cmatrix rho_inv, Carray rho2)
{
    int j,
        s,
        q,
        l,
        i_rho2;

    for (j = 0; j < M; j++)
    {
        for (s = 0; s < M; s++)
        {
            for (q = 0; q < M; q++)
            {
                for (l = 0; l < M; l++)
                {
                    i_rho2 = j + s * M + q * M * M + l * M * M * M;
                    ans += rho_inv[k][j] * rho2[i_rho2] * \
                           conj(Omat[s][n]) * Omat[l][n] * Omat[q][n];
                }
            }
        }
    }
}

void RK4step(MCorbitals psi, MCcoef coef, EqSetup eq, double dt)
{   // Apply 4-th order Runge-Kutta routine given a time step

    long i; // Coeficient Index Counter

    int  k, // Orbital counter
         j, // discretized position counter
         M = psi->Morb,
         Mpos = psi->Mpos,
         Npar = coef->Npar;

    Carray Crhs = carrDef(coef->nc);
    Carray Cnew = carrDef(coef->nc);
    Carray Carg = carrDef(coef->nc);

    Cmatrix Orhs = cmatDef(M, Mpos);
    Cmatrix Onew = cmatDef(M, Mpos);
    Cmatrix Oarg = cmatDef(M, Mpos);

    Cmatrix Ho = cmatdef(M, M);
    Carray Hint = carrDef(M * M * M *M);
    
    /***************               COMPUTE K1               ***************/

    RHSforRK4(Npar, M, coef->NCmat, coef->IF, Mpos, psi->dx,
              eq, coef->C, psi->O, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < coef->nc; i++)
    {   // Add K1 contribution
        Cnew[i] = Crhs[i];
        // Prepare next argument to compute K2
        Carg[i] = coef->C[i] + Crhs[i] * 0.5 * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < N; j++)
        {   // Add K1 contribution
            Onew[k][j] = Orhs[k][j];
            // Prepare next argument to compute K2
            Oarg[k][j] = psi->O[k][j] + Orhs[k][j] * 0.5 * dt;
        }
    }

    /***************               COMPUTE K2               ***************/

    RHSforRK4(Npar, M, coef->NCmat, coef->IF, Mpos, psi->dx,
              eq, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < coef->nc; i++)
    {   // Add K2 contribution
        Cnew[i] += 2 * Crhs[i];
        // Prepare next argument to compute K3
        Carg[i] = coef->C[i] + Crhs[i] * 0.5 * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < N; j++)
        {   // Add K2 contribution
            Onew[k][j] += 2 * Orhs[k][j];
            // Prepare next argument to compute K3
            Oarg[k][j] = psi->O[k][j] + Orhs[k][j] * 0.5 * dt;
        }
    }

    /***************               COMPUTE K3               ***************/

    RHSforRK4(Npar, M, coef->NCmat, coef->IF, Mpos, psi->dx,
              eq, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < coef->nc; i++)
    {   // Add K3 contribution
        Cnew[i] += 2 * Crhs[i];
        // Prepare next argument to compute K4
        Carg->C[i] = coef->C[i] + Crhs[i] * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < N; j++)
        {   // Add K3 contribution
            Onew[k][j] += 2 * Orhs[k][j];
            // Prepare next argument to compute K4
            Oarg->O[k][j] = psi->O[k][j] + Orhs[k][j] * dt;
        }
    }

    /***************               COMPUTE K4               ***************/

    RHSforRK4(Npar, M, coef->NCmat, coef->IF, Mpos, psi->dx,
              eq, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < coef->nc; i++)
    {   // Add K4 contribution
        Cnew[i] += Crhs[i];
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < N; j++)
        {   // Add K4 contribution
            Onew[k][j] += Orhs[k][j];
        }
    }

    // Until now ?new  holds the sum K1 + 2 * K2 + 2 * K3 + K4
    // from the Fourth order Runge-Kutta algorithm. Therefore:
    
    for (i = 0; i < coef->nc; i++) coef->C[i] = coef->C[i] + Cnew[i] * dt / 6;

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < N; j++)
        {   // Add K4 contribution
            psi->O[k][j] += psi->O[k][j] + Onew[k][j] * dt / 6;
        }
    }

    free(Cnew);
    free(Crhs);
    free(Carg);
    
    cmatFree(M, Orhs);
    cmatFree(M, Onew);
    cmatFree(M, Oarg);
}

void RHSforRK4(int N, int M, long ** NCmat, int ** IF,
               int Mpos, double dx, EqSetup eq,
               Carray C, Cmatrix Orb,
               Cmatrix Ho, Carray Hint,
               Carray newC, Cmatrix newOrb)
{
    long i, // index of coeficient
         j,
         nc = NCmat[N][M]; // Total # of coeficients

    int  k, // enumerate orbitals
         l,
         s,
         q;

    /* ================================================================== */
    /*                                                                    */
    /*                        Auxiliary variables                         */
    /*                                                                    */
    /* ================================================================== */

    int  M2 = M * M,
         M3 = M * M * M; // To assist accessing two-body matrix elements

    int * v;             // Occupation vector on each iteration

    double sqrtOf;       // Factor from creation/annihilation operator

    double complex rhsI, // line step value of right-hand-side for C
                   Proj; // Hold value of projector contribution for orbitals

    rhsI = 0;
    Proj = 0;

    /* ================================================================== */
    /*                                                                    */
    /*                       Setup Density Matrices                       */
    /*                                                                    */
    /* ================================================================== */

    Cmatrix rho = cmatDef(M, M);
    Cmatrix rho_inv = cmatDef(M, M);
    Carray rho2 = carrDef(M * M * M * M);

    OBrho(N, M, NCmat, IF, C, rho);
    TBrho(N, M, NCmat, IF, C, rho2);

    HermitianInv(M, rho, rho_inv); // Needed in nonlinear orbital part

    /* ================================================================== */
    /*                                                                    */
    /*                Setup Ho and Hint given the orbitals                */
    /*                                                                    */
    /* ================================================================== */

    SetupHo(M, Mpos, Orb, dx, eq->a2, eq->a1, eq->V, Ho);
    SetupHint(M, Mpos, Orb, dx, eq->inter, Hint);

    /* ================================================================== */
    /*                                                                    */
    /*                      Apply RHS of coeficients                      */
    /*                                                                    */
    /* ================================================================== */

    #pragma omp parallel firstprivate(N, M, M2, M3) \
    private(i, j, k, l, s, q, rhsI, sqrtOf, v)
    {

    v = (int * ) malloc(M * sizeof(int));

    #pragma omp for
    for (i = 0; i < nc; i++)
    {
        rhsI = 0;

        for (k = 0; k < M; k++) v[k] = IF[i][k];
    
        /* ============================================================== */
        /*                                                                */
        /*                      One-body contribution                     */
        /*                                                                */
        /* ============================================================== */

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

        /* ============================================================== */
        /*                                                                */
        /*                      Two-body contribution                     */
        /*                                                                */
        /* ============================================================== */

        /****************************** Rule 1 ******************************/

        for (k = 0; k < M; k++)
        {
            sqrtOf = v[k] * (v[k] - 1)
            rhsI += Hint[k + M * k + M2 * k + M3 * k] * C[i] * sqrtOf;
        }
        
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
                    for (l = 0; l < M; l ++){
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
                    }   // Finish l
                }       // Finish q
            }           // Finish s
        }               // Finish k

        newC[i] = - I * rhsI;
    }
    
    free(v);

    } // End of parallel region

    /* ================================================================== */
    /*                                                                    */
    /*                        Apply RHS of Orbitals                       */
    /*                                                                    */
    /* ================================================================== */

    for (k = 0; k < M; k++)
    {   // Take k orbital
        Proj = 0;
        for (j = 0; j < Mpos; j++)
        {   // At discretized position j
            for (s = 0; s < M; s++)
            {   // projection over 's' orbitals
                Proj += Ho[s][k] * Orb[s][j]; 
                Proj += Proj_Hint(M, k, s, rho_inv, rho2, Hint) * Orb[s][j];
            }
            newOrb[k][j] = -I * (Proj + \
                           eq->inter * NonLinear(M, k, j, Orb, rho_inv, rho2));
        }
    }

    free(rho2);
    cmatFree(M, rho);
    cmatFree(M, rho_inv);

    /**********                  END OF ROUTINE                  **********/
}
