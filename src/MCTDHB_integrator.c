#include "../include/MCTDHB_integrator.h"

MCTDHBsetup AllocMCTDHBdata (
        int Npar, 
        int Morb,
        int Mpos,
        double xi,
        double xf,
        double a2,
        double inter,
        double * V,
        double complex a1   )
{   // Configure and return the pointer to MCTDHB structure
    MCTDHBsetup MC = (MCTDHBsetup) malloc(sizeof(struct _MCTDHBsetup));
    MC->Npar = Npar;
    MC->Morb = Morb;
    MC->Mpos = Mpos;
    MC->xi = xi;
    MC->xf = xf;
    MC->dx = (xf - xi) / (Mpos - 1);
    MC->inter = inter;
    MC->a2 = a2;
    MC->a1 = a1;
    MC->V = V;
    MC->nc = NC(Npar, Morb);
    MC->NCmat = MountNCmat(Npar, Morb);
    MC->IF = MountFocks(Npar, Morb, MC->NCmat);
    return MC;
}

void EraseMCTDHBdata (MCTDHBsetup MC)
{
    long i;
    for (i = 0; i < MC->nc; i++) free(MC->IF[i]);
    free(MC->IF);
    for (i = 0; i <= MC->Npar; i++) free(MC->NCmat[i]);
    free(MC->NCmat);
    free(MC->V);
    free(MC);
}

void SetupHo (int Morb, int Mpos, Cmatrix Omat, double dx, double a2,
              double complex a1, Rarray V, Cmatrix Ho)
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
        for (j = 0; j < Morb; j++)
        {
            dxCyclic(Mpos, Omat[j], dx, ddxj);
            for (k = 0; k < Mpos; k++)
            {
                toInt[k] = - a2 * conj(ddxi[k]) * ddxj[k]    \
                           + a1 * conj(Omat[i][k]) * ddxj[k] \
                           + V[k] * conj(Omat[i][k]) * Omat[j][k];
            }
            Ho[i][j] = Csimps(Mpos, toInt, dx);
            // Ho[j][i] = conj(Ho[i][j]);
        }
    }

    free(ddxi); free(ddxj); free(toInt);
}

void SetupHint (int Morb, int Mpos, Cmatrix Omat, double dx,
                double inter, Carray Hint)
{
    int i,
        k,
        s,
        q,
        l,
        M = Morb,
        M2 = Morb * Morb,
        M3 = Morb * Morb * Morb;

    double complex Integral;

    Carray toInt = carrDef(Mpos);

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
                                   Omat[q][i] * Omat[l][i];
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

double complex Proj_Hint(int M, int k, int i, Cmatrix rho_inv, Carray rho2,
                         Carray Hint)
{   // k and i enumerate orbitals (see lecture notes)
    int j,
        s,
        q,
        l,
        i_rho2,
        i_Hint;

    double complex ans = 0;

    for (j = 0; j < M; j++)
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

double complex NonLinear(int M, int k, int n, Cmatrix Omat, Cmatrix rho_inv,
                         Carray rho2)
{   // k enumerate orbital and n a discretized position
    int j,
        s,
        q,
        l,
        i_rho2;

    double complex ans = 0;

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

    return ans;
}

void RK4step(MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt)
{   // Apply 4-th order Runge-Kutta routine given a time step

    long i; // Coeficient Index Counter

    int  k, // Orbital counter
         j, // discretized position counter
         M = MC->Morb,
         Mpos = MC->Mpos,
         Npar = MC->Npar;

    Carray Crhs = carrDef(MC->nc);
    Carray Cnew = carrDef(MC->nc);
    Carray Carg = carrDef(MC->nc);

    Cmatrix Orhs = cmatDef(M, Mpos);
    Cmatrix Onew = cmatDef(M, Mpos);
    Cmatrix Oarg = cmatDef(M, Mpos);

    Cmatrix  Ho = cmatDef(M, M);
    Carray Hint = carrDef(M * M * M *M);



    /* __________________________________________________________________ *
     *                             COMPUTE K1                             *
     *                          ----------------                          */

    RHSforRK4(MC, C, Orb, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K1 contribution
        Cnew[i] = Crhs[i];
        // Prepare next argument to compute K2
        Carg[i] = C[i] + Crhs[i] * 0.5 * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K1 contribution
            Onew[k][j] = Orhs[k][j];
            // Prepare next argument to compute K2
            Oarg[k][j] = Orb[k][j] + Orhs[k][j] * 0.5 * dt;
        }
    }


    /* __________________________________________________________________ *
     *                             COMPUTE K2                             *
     *                          ----------------                          */

    RHSforRK4(MC, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K2 contribution
        Cnew[i] += 2 * Crhs[i];
        // Prepare next argument to compute K3
        Carg[i] = C[i] + Crhs[i] * 0.5 * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K2 contribution
            Onew[k][j] += 2 * Orhs[k][j];
            // Prepare next argument to compute K3
            Oarg[k][j] = Orb[k][j] + Orhs[k][j] * 0.5 * dt;
        }
    }




    /* __________________________________________________________________ *
     *                             COMPUTE K3                             *
     *                          ----------------                          */

    RHSforRK4(MC, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K3 contribution
        Cnew[i] += 2 * Crhs[i];
        // Prepare next argument to compute K4
        Carg[i] = C[i] + Crhs[i] * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K3 contribution
            Onew[k][j] += 2 * Orhs[k][j];
            // Prepare next argument to compute K4
            Oarg[k][j] = Orb[k][j] + Orhs[k][j] * dt;
        }
    }



    /* __________________________________________________________________ *
     *                             COMPUTE K3                             *
     *                          ----------------                          */

    RHSforRK4(MC, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K4 contribution
        Cnew[i] += Crhs[i];
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K4 contribution
            Onew[k][j] += Orhs[k][j];
        }
    }

    // Until now ?new  holds the sum K1 + 2 * K2 + 2 * K3 + K4
    // from the Fourth order Runge-Kutta algorithm. Therefore:

    for (i = 0; i < MC->nc; i++)
    {   // Update Coeficients
        C[i] = C[i] + Cnew[i] * dt / 6;
    }

    for (k = 0; k < M; k++)
    {   // Update Orbitals
        for (j = 0; j < Mpos; j++)
        {
            Orb[k][j] = Orb[k][j] + Onew[k][j] * dt / 6;
        }
    }

    free(Cnew);
    free(Crhs);
    free(Carg);
    
    cmatFree(M, Orhs);
    cmatFree(M, Onew);
    cmatFree(M, Oarg);

    cmatFree(M, Ho);
    free(Hint);
}

void IRK4step(MCTDHBsetup MC, Cmatrix Orb, Carray C, double complex dt)
{   // Apply 4-th order Runge-Kutta routine given a time step

    long i; // Coeficient Index Counter

    int  k, // Orbital counter
         j, // discretized position counter
         M = MC->Morb,
         Mpos = MC->Mpos,
         Npar = MC->Npar;

    Carray Crhs = carrDef(MC->nc);
    Carray Cnew = carrDef(MC->nc);
    Carray Carg = carrDef(MC->nc);

    Cmatrix Orhs = cmatDef(M, Mpos);
    Cmatrix Onew = cmatDef(M, Mpos);
    Cmatrix Oarg = cmatDef(M, Mpos);

    Cmatrix  Ho = cmatDef(M, M);
    Carray Hint = carrDef(M * M * M *M);



    /* __________________________________________________________________ *
     *                             COMPUTE K1                             *
     *                          ----------------                          */

    RHSforRK4(MC, C, Orb, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K1 contribution
        Cnew[i] = Crhs[i];
        // Prepare next argument to compute K2
        Carg[i] = C[i] + Crhs[i] * 0.5 * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K1 contribution
            Onew[k][j] = Orhs[k][j];
            // Prepare next argument to compute K2
            Oarg[k][j] = Orb[k][j] + Orhs[k][j] * 0.5 * dt;
        }
    }



    /* __________________________________________________________________ *
     *                             COMPUTE K2                             *
     *                          ----------------                          */

    RHSforRK4(MC, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K2 contribution
        Cnew[i] += 2 * Crhs[i];
        // Prepare next argument to compute K3
        Carg[i] = C[i] + Crhs[i] * 0.5 * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K2 contribution
            Onew[k][j] += 2 * Orhs[k][j];
            // Prepare next argument to compute K3
            Oarg[k][j] = Orb[k][j] + Orhs[k][j] * 0.5 * dt;
        }
    }



    /* __________________________________________________________________ *
     *                             COMPUTE K3                             *
     *                          ----------------                          */

    RHSforRK4(MC, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K3 contribution
        Cnew[i] += 2 * Crhs[i];
        // Prepare next argument to compute K4
        Carg[i] = C[i] + Crhs[i] * dt;
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K3 contribution
            Onew[k][j] += 2 * Orhs[k][j];
            // Prepare next argument to compute K4
            Oarg[k][j] = Orb[k][j] + Orhs[k][j] * dt;
        }
    }



    /* __________________________________________________________________ *
     *                             COMPUTE K3                             *
     *                          ----------------                          */

    RHSforRK4(MC, Carg, Oarg, Ho, Hint, Crhs, Orhs);

    for (i = 0; i < MC->nc; i++)
    {   // Add K4 contribution
        Cnew[i] += Crhs[i];
    }

    for (k = 0; k < M; k++)
    {
        for (j = 0; j < Mpos; j++)
        {   // Add K4 contribution
            Onew[k][j] += Orhs[k][j];
        }
    }

    // Until now ?new  holds the sum K1 + 2 * K2 + 2 * K3 + K4
    // from the Fourth order Runge-Kutta algorithm. Therefore:

    for (i = 0; i < MC->nc; i++)
    {   // Update Coeficients
        C[i] = C[i] + Cnew[i] * dt / 6;
    }

    for (k = 0; k < M; k++)
    {   // Update Orbitals
        for (j = 0; j < Mpos; j++)
        {
            Orb[k][j] += Orb[k][j] + Onew[k][j] * dt / 6;
        }
    }

    free(Cnew);
    free(Crhs);
    free(Carg);
    
    cmatFree(M, Orhs);
    cmatFree(M, Onew);
    cmatFree(M, Oarg);

    cmatFree(M, Ho);
    free(Hint);
}

void RHSforRK4(MCTDHBsetup MC, Carray C, Cmatrix Orb,
               Cmatrix Ho, Carray Hint,
               Carray newC, Cmatrix newOrb)
{
    long // Index of coeficients
        i,
        j,
        nc = MC->nc;

    int // enumerate orbitals
        k,
        l,
        s,
        q;

    /* ================================================================== */
    /*                                                                    */
    /*                        Auxiliary variables                         */
    /*                                                                    */
    /* ================================================================== */

    int  M  = MC->Morb,
         N  = MC->Npar,
         M2 = M * M,
         M3 = M * M * M,
         Mpos  = MC->Mpos,
         ** IF = MC->IF;
    
    long ** NCmat = MC->NCmat;

    double dx = MC->dx;

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

    l = HermitianInv(M, rho, rho_inv); // Needed in nonlinear orbital part



    /* ================================================================== */
    /*                                                                    */
    /*                Setup Ho and Hint given the orbitals                */
    /*                                                                    */
    /* ================================================================== */

    SetupHo(M, Mpos, Orb, dx, MC->a2, MC->a1, MC->V, Ho);
    SetupHint(M, Mpos, Orb, dx, MC->inter, Hint);



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
            sqrtOf = v[k] * (v[k] - 1);
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
        for (j = 0; j < Mpos; j++)
        {   // At discretized position j
            Proj = 0;
            for (s = 0; s < M; s++)
            {   // projection over 's' orbitals
                Proj += Ho[s][k] * Orb[s][j]; 
                Proj += Proj_Hint(M, k, s, rho_inv, rho2, Hint) * Orb[s][j];
            }
            newOrb[k][j] = -I * (MC->inter * NonLinear(M, k, j, Orb, rho_inv, rho2) - Proj);
        }
    }

    free(rho2);
    cmatFree(M, rho);
    cmatFree(M, rho_inv);

    /**********                  END OF ROUTINE                  **********/
}

void LinearPartSM(MCTDHBsetup MC, CCSmat rhs_mat, Carray upper, Carray lower,
                Carray mid, Cmatrix Orb)
{
    int k,
        size = MC->Mpos - 1;

    Carray rhs = carrDef(size);

    for (k = 0; k < MC->Morb; k++)
    {   // For each orbital k solve a tridiagonal system obtained by CN
        CCSvec(size, rhs_mat->vec, rhs_mat->col, rhs_mat->m, Orb[k], rhs);
        triCyclicSM(size, upper, lower, mid, rhs, Orb[k]);
        Orb[k][size] = Orb[k][0]; // the boundary
    }

    free(rhs);
}

void MCTDHB_time_evolution(MCTDHBsetup MC, Cmatrix Orb, Carray C, double dt,
        int Nsteps, int cyclic)
{
    int i,
        Mpos = MC->Mpos;

    double dx = MC->dx,
           a2 = MC->a2;

    double complex a1 = MC->a1;

    // used to store matrix elements of linear part
    Carray upper = carrDef(Mpos - 1);
    Carray lower = carrDef(Mpos - 1);
    Carray mid   = carrDef(Mpos - 1);



    /*                                                               *
     *      Setup Right-Hand-Side matrix of linear part of PDE       *
     *      --------------------------------------------------       */



    // fill main diagonal (use upper as auxiliar array)
    carrFill(Mpos - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(Mpos - 1, upper, dt, MC->V, mid);

    // fill upper diagonal
    carrFill(Mpos - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[Mpos-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[Mpos-2] = 0;                                        }

    // fill lower diagonal
    carrFill(Mpos - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[Mpos-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[Mpos-2] = 0;                                        }

    // Store in CCS format the RHS of discretized system of equations
    CCSmat rhs_mat = CyclicToCCS(Mpos - 1, upper, lower, mid);



    /*                                                               *
     *                Setup Cyclic tridiagonal matrix                *
     *                -------------------------------                */



    // fill main diagonal (use upper as auxiliar array)
    carrFill(Mpos - 1, a2 * dt / dx / dx + I, upper);
    rcarrUpdate(Mpos - 1, upper, -dt, MC->V, mid);

    // fill upper diagonal
    carrFill(Mpos - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[Mpos-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[Mpos-2] = 0;                                          }

    // fill lower diagonal
    carrFill(Mpos - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[Mpos-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[Mpos-2] = 0;                                          }

    for (i = 0; i < Nsteps; i++)
    {
        RK4step(MC, Orb, C, dt / 2);
        LinearPartSM(MC, rhs_mat, upper, lower, mid, Orb);
        RK4step(MC, Orb, C, dt / 2);
    }

    CCSFree(rhs_mat);
    free(upper);
    free(lower);
    free(mid);
}

void MCTDHB_itime_evolution(MCTDHBsetup MC, Cmatrix Orb, Carray C, double dT,
        int Nsteps, int cyclic)
{
    int i,
        k,
        Mpos = MC->Mpos;

    double dx = MC->dx,
           a2 = MC->a2;

    double complex a1 = MC->a1,
                   dt = - I * dT;

    // used to store matrix elements of linear part
    Carray upper = carrDef(Mpos - 1);
    Carray lower = carrDef(Mpos - 1);
    Carray mid   = carrDef(Mpos - 1);



    /*                                                               *
     *      Setup Right-Hand-Side matrix of linear part of PDE       *
     *      --------------------------------------------------       */



    // fill main diagonal (use upper as auxiliar array)
    carrFill(Mpos - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(Mpos - 1, upper, dt, MC->V, mid);

    // fill upper diagonal
    carrFill(Mpos - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[Mpos-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[Mpos-2] = 0;                                        }

    // fill lower diagonal
    carrFill(Mpos - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[Mpos-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[Mpos-2] = 0;                                        }

    // Store in CCS format the RHS of discretized system of equations
    CCSmat rhs_mat = CyclicToCCS(Mpos - 1, upper, lower, mid);



    /*                                                               *
     *                Setup Cyclic tridiagonal matrix                *
     *                -------------------------------                */



    // fill main diagonal (use upper as auxiliar array)
    carrFill(Mpos - 1, a2 * dt / dx / dx + I, upper);
    rcarrUpdate(Mpos - 1, upper, -dt, MC->V, mid);

    // fill upper diagonal
    carrFill(Mpos - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[Mpos-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[Mpos-2] = 0;                                          }

    // fill lower diagonal
    carrFill(Mpos - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[Mpos-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[Mpos-2] = 0;                                          }

    for (i = 0; i < Nsteps; i++)
    {
        IRK4step(MC, Orb, C, dt / 2);
        LinearPartSM(MC, rhs_mat, upper, lower, mid, Orb);
        IRK4step(MC, Orb, C, dt / 2);
        // Compared to real time we need additional renormalization
        for (k = 0; k < MC->Morb; k++)
        {   // Renormalize each orbital
            renormalize(Mpos, Orb[k], dx, 1.0);
        }
        renormalizeVector(MC->nc, C, 1.0);
    }

    CCSFree(rhs_mat);
    free(upper);
    free(lower);
    free(mid);
}
