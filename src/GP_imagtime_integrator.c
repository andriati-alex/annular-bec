#include "../include/GP_imagtime_integrator.h"





void IGPFFT(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, Carray S, Carray E)
{

    int
        i,
        j,
        m;

    m = M - 1; // ignore last point that is the boundary

    double
        freq,       // frequencies in Fourier space
        norm,       // initial norm
        NormStep,   // to renormalize on each time-step
        Idt = - dT; // factor to multiply on exponential after split-step

    double complex
        dt = - I  * dT; // imaginary time-step



    MKL_LONG s; // status of called MKL FFT functions

    Rarray abs2 = rarrDef(M);       // abs square of wave function

    Rarray out  = rarrDef(M);       // hold linear and nonlinear potential

    Carray exp_der = carrDef(m);    // exponential of derivative operator

    Carray stepexp = carrDef(M);    // Exponential of potential

    Carray back_fft = carrDef(m);   // go back to position space

    Carray foward_fft = carrDef(m); // go to frequency space



    /* Initialize the norm and energy of initial guess
     * ------------------------------------------------------------------- */
    carrAbs2(M, S, abs2);
    norm = sqrt(Rsimps(M, abs2, dx));
    E[0] = Functional(M, dx, a2, a1, inter / 2, V, S);
    /* ------------------------------------------------------------------- */



    /* setup descriptor (MKL implementation of FFT)
     * ------------------------------------------------------------------- */
    DFTI_DESCRIPTOR_HANDLE desc;
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, m);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(m));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(m));
    s = DftiCommitDescriptor(desc);
    /* ------------------------------------------------------------------- */



    /* Fourier Frequencies to do the exponential of derivative operator
     * ------------------------------------------------------------------- */
    for (i = 0; i < m; i++)
    {
        if (i <= (m - 1) / 2) { freq = (2 * PI * i) / (m * dx);       }
        else                  { freq = (2 * PI * (i - m)) / (m * dx); }
        // exponential of derivative operators
        exp_der[i] = cexp((-1) * Idt * a2 * freq * freq);
    }
    /* ------------------------------------------------------------------- */



    /*   Apply Split step and solve separately nonlinear and linear part   */
    /*   ===============================================================   */

    for (i = 0; i < N; i++) {
        // Apply exponential ofone-body potential and nonlinear part
        rarrUpdate(M, V, inter, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);
        carrMultiply(m, stepexp, S, foward_fft);

        s = DftiComputeForward(desc, foward_fft);       // go to momentum space
        carrMultiply(m, exp_der, foward_fft, back_fft); // apply derivatives
        s = DftiComputeBackward(desc, back_fft);        // back to real space
        carrMultiply(m, stepexp, back_fft, S);
        S[m] = S[0];

        carrAbs2(M, S, abs2);

        // Renormalization
        NormStep = norm / sqrt(Rsimps(M, abs2, dx));
        for (j = 0; j < M; j++) S[j] = NormStep * S[j];

        // Energy
        E[i + 1] = Functional(M, dx, a2, a1, inter / 2, V, S);
    }

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(stepexp);
    free(foward_fft);
    free(back_fft);
    free(abs2);
    free(out);
}





void IGPCNSM(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, Carray E)
{
    unsigned int
        i,
        j;

    double
        norm,       // initial norm
        NormStep,   // to renormalize at each time-step
        Idt = - dT; // factor that multiplies in split-step exponentials

    double complex
        dt = - I  * dT; // pure imaginary time-step

    // Used to apply nonlinear part
    Carray stepexp = carrDef(M);
    Carray linpart = carrDef(M);
    Rarray abs2    = rarrDef(M); // abs square of wave function

    // used to store matrix elements of linear part
    Carray upper = carrDef(M - 1);
    Carray lower = carrDef(M - 1);
    Carray mid   = carrDef(M - 1);
    Carray rhs   = carrDef(M - 1); // RHS of linear system from CN

    carrAbs2(M, S, abs2);
    norm = sqrt(Rsimps(M, abs2, dx));
    E[0] = Functional(M, dx, a2, a1, inter / 2, V, S);



    /*                 ****************************                 */
    /*                 Setup Right-Hand-Side matrix                 */
    /*                 ****************************                 */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M - 1, upper, dt, V, mid);

    // fill upper diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                        }

    // fill lower diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                        }

    // Store in CCS format
    CCSmat rhs_mat = CyclicToCCS(M - 1, upper, lower, mid);



    /*                *******************************                */
    /*                Setup Cyclic tridiagonal matrix                */
    /*                *******************************                */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(M - 1, upper, -dt, V, mid);

    // fill upper diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                          }

    // fill lower diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                          }



    /* Apply Split step and solve separately nonlinear and linear part */
    /* *************************************************************** */



    for (i = 0; i < N; i++) {
        // Apply exponential with nonlinear part
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        // Update solution (solution of linear part in stepexp)
        carrMultiply(M, linpart, stepexp, S);

        carrAbs2(M, S, abs2);

        // Renormalize
        NormStep = norm / sqrt(Rsimps(M, abs2, dx));
        for (j = 0; j < M; j++) S[j] = NormStep * S[j];
        
        // Energy
        E[i + 1] = Functional(M, dx, a2, a1, inter / 2, V, S);
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





void IGPCNLU(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, Carray E)
{
    unsigned int
        i,
        j;

    double
        norm,
        NormStep,   // to renormalize at each time-step
        Idt = - dT; // factor that multiplies in split-step exponentials

    double complex
        dt = - I  * dT;

    // Used to apply nonlinear part
    Carray stepexp = carrDef(M);
    Carray linpart = carrDef(M);
    Rarray abs2    = rarrDef(M); // abs square of wave function

    // used to store matrix elements of linear part
    Carray upper = carrDef(M - 1);
    Carray lower = carrDef(M - 1);
    Carray mid   = carrDef(M - 1);
    Carray rhs   = carrDef(M - 1);

    carrAbs2(M, S, abs2);
    norm = sqrt(Rsimps(M, abs2, dx));
    E[0] = Functional(M, dx, a2, a1, inter / 2, V, S);



    /*                 ****************************                 */
    /*                 Setup Right-Hand-Side matrix                 */
    /*                 ****************************                 */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M - 1, upper, dt, V, mid);

    // fill upper diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                        }

    // fill lower diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                        }

    // Store in CCS format
    CCSmat rhs_mat = CyclicToCCS(M - 1, upper, lower, mid);



    /*                *******************************                */
    /*                Setup Cyclic tridiagonal matrix                */
    /*                *******************************                */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(M - 1, upper, -dt, V, mid);

    // fill upper diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                          }

    // fill lower diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                          }



    /* Apply Split step and solve separately nonlinear and linear part */
    /* *************************************************************** */



    for (i = 0; i < N; i++) {
        // Apply exponential with nonlinear part
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicLU(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        // Update solution (solution of linear part in stepexp)
        carrMultiply(M, linpart, stepexp, S);

        carrAbs2(M, S, abs2);

        // Renormalize
        NormStep = norm / sqrt(Rsimps(M, abs2, dx));
        for (j = 0; j < M; j++) S[j] = NormStep * S[j];
        
        // Energy
        E[i + 1] = Functional(M, dx, a2, a1, inter / 2, V, S);
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



void NonLinearIDDT(int M, double t, Carray Psi, Carray inter, Carray Dpsi)
{   // The Right-Hand-Side of derivative of non-linear part
    // As a extra argument take the interaction strength

    int i;

    // sum real and imaginary part squared to get the modulus squared
    double mod_squared;

    for (i = 0; i < M; i++)
    {
        mod_squared  = creal(Psi[i]) * creal(Psi[i]);
        mod_squared += cimag(Psi[i]) * cimag(Psi[i]);
        Dpsi[i] = - inter[0] * mod_squared * Psi[i];
    }
}



void IGPCNSMRK4(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, Carray E)
{   // Evolve Gross-Pitaevskii using 4-th order Runge-Kutta
    // to deal with nonlinear part.
    
    int i,
        j;

    double
        norm,
        NormStep;

    double complex
        dt,
        interv[1];

    dt = - I * dT; // Set time to pure imaginary

    interv[0] = inter;
    
    Rarray abs2 = rarrDef(M); // abs square of wave function

    // Used to apply nonlinear part
    Carray linpart = carrDef(M);

    // used to store matrix elements of linear part
    Carray upper = carrDef(M - 1);
    Carray lower = carrDef(M - 1);
    Carray mid   = carrDef(M - 1);
    Carray rhs   = carrDef(M - 1);
    
    // Initial Energy and norm
    carrAbs2(M, S, abs2);
    norm = sqrt(Rsimps(M, abs2, dx));
    E[0] = Functional(M, dx, a2, a1, inter / 2, V, S);



    /*                 ****************************                 */
    /*                 Setup Right-Hand-Side matrix                 */
    /*                 ****************************                 */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M - 1, upper, dt, V, mid);

    // fill upper diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                        }

    // fill lower diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                        }

    // Store in CCS format
    CCSmat rhs_mat = CyclicToCCS(M - 1, upper, lower, mid);



    /*                *******************************                */
    /*                Setup Cyclic tridiagonal matrix                */
    /*                *******************************                */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(M - 1, upper, -dt, V, mid);

    // fill upper diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                          }

    // fill lower diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                          }



    /* Apply Split step and solve separately nonlinear and linear part */
    /* *************************************************************** */



    for (i = 0; i < N; i++)
    {
        RK4step(M, dT/2, 0, S, interv, linpart, NonLinearIDDT);
        
        // Solve linear part (nabla ^ 2 part)
        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        RK4step(M, dT/2, 0, linpart, interv, S, NonLinearIDDT);

        carrAbs2(M, S, abs2);
        
        // Renormalize
        NormStep = norm / sqrt(Rsimps(M, abs2, dx));
        for (j = 0; j < M; j++) S[j] = NormStep * S[j];
        
        // Energy
        E[i + 1] = Functional(M, dx, a2, a1, inter / 2, V, S);
    }

    free(linpart);
    free(upper);
    free(lower);
    free(abs2);
    free(mid);
    free(rhs);
    CCSFree(rhs_mat);
}
