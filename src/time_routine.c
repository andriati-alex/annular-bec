#include "../include/time_routine.h"

void CNsm(int M, int N, double dx, double dt, double a2, double complex a1,
          double inter, Rarray V, int cyclic, Cmatrix S)
{
    unsigned int i, j;
    double complex Idt = 0.0 - dt * I;

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
    carrFill(M, a2 * dt / dx /dx + I, upper);
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
        carrAbs2(M, S[i], abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S[i], linpart);

        // Solve linear part
        CCSvec(M, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M, upper, lower, mid, rhs, linpart);

        // Apply again the exponential of nonlinear part
        carrMultiply(M, linpart, stepexp, S[i+1]);

        /* IMPROVEMENT
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the nonlinear part.
         *
         */

        // Trapezium rule in the nonlinear exponential part
        for (j = 0; j < M; j++) {
            abs2[j] += creal(S[i+1][j]) * creal(S[i+1][j]);
            abs2[j] += cimag(S[i+1][j]) * cimag(S[i+1][j]);
        }

        // Extra factor 1/2 due to trapezium rule
        rcarrExp(M, inter * Idt / 4, abs2, stepexp);
        carrMultiply(M, stepexp, S[i], linpart);

        CCSvec(M, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M, upper, lower, mid, rhs, linpart);

        carrMultiply(M, linpart, stepexp, S[i+1]);
    }

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    CCSFree(rhs_mat);
}



void CNlu(int M, int N, double dx, double dt, double a2, double complex a1,
          double inter, Rarray V, int cyclic, Cmatrix S)
{
    unsigned int i, j;
    double complex Idt = 0.0 - dt * I;

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
    carrFill(M, a2 * dt / dx /dx + I, upper);
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
        carrAbs2(M, S[i], abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S[i], linpart);

        // Solve linear part
        CCSvec(M, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicLU(M, upper, lower, mid, rhs, linpart);

        // Update solution (solution of linear part in stepexp)
        carrMultiply(M, linpart, stepexp, S[i+1]);
        
        /* IMPROVEMENT
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the nonlinear part.
         *
         */

        // Trapezium rule in the nonlinear exponential part
        for (j = 0; j < M; j++) {
            abs2[j] += creal(S[i+1][j]) * creal(S[i+1][j]);
            abs2[j] += cimag(S[i+1][j]) * cimag(S[i+1][j]);
        }

        // Extra factor 1/2 due to trapezium rule
        rcarrExp(M, inter * Idt / 4, abs2, stepexp);
        carrMultiply(M, stepexp, S[i], linpart);

        CCSvec(M, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicLU(M, upper, lower, mid, rhs, linpart);

        carrMultiply(M, linpart, stepexp, S[i+1]);
    }

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    CCSFree(rhs_mat);
}



void spectral(int M, int N, double dx, double dt, double a2, double complex a1,
              double inter, Rarray V, Cmatrix S)
{
    int i, j;
    double freq;
    double complex Idt = 0.0 - dt * I;

    // mkl_set_num_threads(mkl_get_max_threads() / 2);

    MKL_LONG s; // status of called MKL FFT functions

    Rarray abs2 = rarrDef(M); // abs square of wave function
    Rarray out  = rarrDef(M); // output of linear and nonlinear potential

    Carray exp_der = carrDef(M); // exponential of derivative operator

    Carray stepexp = carrDef(M);     // Exponential of potential
    Carray back_fft = carrDef(M);    // back to position space

    /***************** setup descriptor (MKL implementation) ****************/

    DFTI_DESCRIPTOR_HANDLE desc;
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, M);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(M));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(M));
    s = DftiCommitDescriptor(desc);

    /************************************************************************/

    /* setup Fourier Frequencies and the exponential of derivative operator */
    for (i = 0; i < M; i++) {
        if (i <= (M - 1) / 2) { freq = (2 * PI * i) / (M * dx);       }
        else                  { freq = (2 * PI * (i - M)) / (M * dx); }
        // exponential of derivative operators
        exp_der[i] = cexp((-1) * Idt * a2 * freq * freq);
    }
    
    /* Apply Split step and solve separately nonlinear and linear part */
    /*******************************************************************/

    for (i = 0; i < N; i++) {
        // Apply exponential of potential together nonlinear part
        carrAbs2(M, S[i], abs2);
        rarrUpdate(M, V, inter, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);
        carrMultiply(M, stepexp, S[i], S[i+1]);

        s = DftiComputeForward(desc, S[i+1]);       // go to momentum space
        carrMultiply(M, exp_der, S[i+1], back_fft); // apply derivatives
        s = DftiComputeBackward(desc, back_fft);    // back to real space
        carrMultiply(M, stepexp, back_fft, S[i+1]);
        
        /* IMPROVEMENT
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the nonlinear part.
         *
         */

        // Trapezium rule in the nonlinear exponential part
        for (j = 0; j < M; j++) {
            abs2[j] += creal(S[i+1][j]) * creal(S[i+1][j]);
            abs2[j] += cimag(S[i+1][j]) * cimag(S[i+1][j]);
        }
        
        rarrUpdate(M, V, inter / 2, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);
        carrMultiply(M, stepexp, S[i], S[i+1]);

        s = DftiComputeForward(desc, S[i+1]);       // go to momentum space
        carrMultiply(M, exp_der, S[i+1], back_fft); // apply derivatives
        s = DftiComputeBackward(desc, back_fft);    // back to real space
        carrMultiply(M, stepexp, back_fft, S[i+1]);
    }

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(stepexp);
    free(back_fft);
    free(abs2);
    free(out);
}



void lastStep(int M, int N, double dx, double dt, double a2, double complex a1,
              double inter, Rarray V, int cyclic, Carray S)
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
