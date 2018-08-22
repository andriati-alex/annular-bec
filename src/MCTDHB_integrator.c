#include "../include/MCTDHB_integrator.h"

void Coef_HalfStep(int N, int M, int ** IF, long ** NCmat, Cmatrix Ho,
                   Carray Hint, Carray C, double dt, Carray Cnext)
{   // Evolve half time step using Euler-Cauchy method(explicit trapezium)

    long i; // Counter

    Carray Caux = carrDef(NCmat[N][M]);

    RHSofODES(N, M, IF, NCmat, Ho, Hint, C, Caux);

    for (i = 0; i < NCmat[N][M]; i++)
    {   // Caux will be the vector to call the foward step RHS
        Caux[i] = (0.5 * dt) * Caux[i] + C[i];
    }

    // Cnext hold the foward (half)time step RHS of the system, for now
    RHSofODES(N, M, IF, NCmat, Ho, Hint, Caux, Cnext);

    for (i = 0; i < NCmat[N][M]; i++)
    {   // Add up the terms to get the solution at next (half)time step
        Caux[i] = (Caux[i] - C[i]); // The first term we must sum
        Cnext[i] = C[i] + Caux[i] * 0.5 + Cnext[i] * (0.25 * dt);
    }

    free(Caux);
}

void Coef_Step(int N, int M, int ** IF, long ** NCmat, Cmatrix Ho,
                   Carray Hint, Carray C, double dt, Carray Cnext)
{   // Evolve half time step using Euler-Cauchy method(explicit trapezium)

    long i; // Counter

    Carray Karg = carrDef(NCmat[N][M]);
    Carray Kans = carrDef(NCmat[N][M]);

    /***************               COMPUTE K1               ***************/

    RHSofODES(N, M, IF, NCmat, Ho, Hint, C, Kans);

    for (i = 0; i < NCmat[N][M]; i++)
    {   // Add K1 contribution
        Cnext[i] = Kans[i];
        // Prepare next argument to compute K2
        Karg[i] = C[i] + Kans[i] * 0.5 * dt;
    }

    /* ================================================================== */
    /*                                                                    */
    /*                     Need to update Ho and Hint                     */
    /*                                                                    */
    /* ================================================================== */
    
    /***************               COMPUTE K2               ***************/

    RHSofODES(N, M, IF, NCmat, Ho, Hint, Karg, Kans);

    for (i = 0; i < NCmat[N][M]; i++)
    {   // Add K2 contribution
        Cnext[i] += 2 * Kans[i];
        // Prepare next argument to compute K3
        Karg[i] = C[i] + Kans[i] * 0.5 * dt;
    }

    /***************               COMPUTE K3               ***************/

    RHSofODES(N, M, IF, NCmat, Ho, Hint, Karg, Kans);
    
    for (i = 0; i < NCmat[N][M]; i++)
    {   // Add K3 contribution
        Cnext[i] += 2 * Kans[i];
        // Prepare next argument to compute K4
        Karg[i] = C[i] + Kans[i] * dt;
    }
    
    /* ================================================================== */
    /*                                                                    */
    /*                     Need to update Ho and Hint                     */
    /*                                                                    */
    /* ================================================================== */

    /***************               COMPUTE K4               ***************/

    RHSofODES(N, M, IF, NCmat, Ho, Hint, Karg, Kans);

    for (i = 0; i < NCmat[N][M]; i++)
    {   // Add K4 contribution
        Cnext[i] += Kans[i];
    }

    // Until now Cnext holds the sum K1 + 2 * K2 + 2 * K3 + K4
    // from the Fourth order Runge-Kutta algorithm. Therefore:
    
    for (i = 0; i < NCmat[N][M]; i++) Cnext[i] = C[i] + Cnext * dt / 6;

}

void FieldStep(Carray C, Cmatrix rho, Carray rho2)
{
    unsigned int i, j;
    double complex Idt = 0.0 - dt * I;

    Carray Sstep = carrDef(M);

    // Used to apply nonlinear part
    Carray stepexp = carrDef(M);
    Carray linpart = carrDef(M);
    Rarray abs2    = rarrDef(M); // abs square of wave function

    // used to store matrix elements of linear part
    Carray upper = carrDef(M);
    Carray lower = carrDef(M);
    Carray mid   = carrDef(M);
    Carray rhs   = carrDef(M);

    /*                 ****************************                 */
    /*                 Setup Right-Hand-Side matrix                 */
    /*                 ****************************                 */

    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M, upper, dt, V, mid);

    // fill upper diagonal
    carrFill(M, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-1] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-1] = 0;                                        }

    // fill lower diagonal
    carrFill(M, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-1] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-1] = 0;                                        }

    // Store in CCS format
    CCSmat rhs_mat = CyclicToCCS(M, upper, lower, mid);

    /*                *******************************                */
    /*                Setup Cyclic tridiagonal matrix                */
    /*                *******************************                */

    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M, a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M, upper, -dt, V, mid);

    // fill upper diagonal
    carrFill(M, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-1] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-1] = 0;                                          }

    // fill lower diagonal
    carrFill(M, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-1] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-1] = 0;                                          }

    /* Apply Split step and solve separately nonlinear and linear part */
    /* *************************************************************** */

    for (i = 0; i < N; i++) {
        // Apply exponential with nonlinear part
        carrAbs2(M, S, abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M, upper, lower, mid, rhs, linpart);

        // Update solution (solution of linear part in stepexp)
        carrMultiply(M, linpart, stepexp, Sstep);
        
        /* IMPROVEMENT
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the nonlinear part.
         *
         */

        // Trapezium rule in the nonlinear exponential part
        for (j = 0; j < M; j++) {
            abs2[j] += creal(Sstep[j]) * creal(Sstep[j]);
            abs2[j] += cimag(Sstep[j]) * cimag(Sstep[j]);
        }

        // Extra factor 1/2 due to trapezium rule
        rcarrExp(M, inter * Idt / 4, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        CCSvec(M, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M, upper, lower, mid, rhs, linpart);

        carrMultiply(M, linpart, stepexp, S);
    }

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    free(Sstep);
    CCSFree(rhs_mat);
}
