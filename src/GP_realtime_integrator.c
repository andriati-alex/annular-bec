#include "../include/GP_realtime_integrator.h"

void record_step(FILE * f, int M, Carray v)
{
    int j;

    for (j = 0; j < M; j ++)
    {
        if (cimag(v[j]) > 0)
        {
            fprintf(f, "(%.15E+%.15Ej) ", creal(v[j]), cimag(v[j]));
        }
        else
        {
            if (cimag(v[j]) == 0)
                fprintf(f, "(%.15E+%.15Ej) ", creal(v[j]), 0.0);
            else
                fprintf(f, "(%.15E%.15Ej) ", creal(v[j]), cimag(v[j]));
        }
    }

    fprintf(f, "\n");
}




void GPFFT(int M, int N, double dx, double dt, double a2,
     double complex a1, double inter, Rarray V, Carray S)
{
    int
        i,
        j;

    double
        freq;

    double complex
        Idt = 0.0 - dt * I;



    MKL_LONG s; // status of called MKL FFT functions

    Rarray abs2 = rarrDef(M);    // abs square of wave function

    Rarray out  = rarrDef(M);    // output of linear and nonlinear potential

    Carray exp_der = carrDef(M); // exponential of derivative operator

    Carray stepexp = carrDef(M);    // Exponential of potential

    Carray foward_fft = carrDef(M); // go to frequency space

    Carray back_fft = carrDef(M);   // back to position space
    
    Carray Sstep = carrDef(M);      // hold one step to integrate by trapezium



    /* setup descriptor (MKL implementation of FFT)
     * -------------------------------------------------------------------- */
    DFTI_DESCRIPTOR_HANDLE desc;
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, M);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(M));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(M));
    s = DftiCommitDescriptor(desc);
    /* -------------------------------------------------------------------- */



    /* setup Fourier Frequencies and the exponential of derivative operator
     * -------------------------------------------------------------------- */
    for (i = 0; i < M; i++) {
        if (i <= (M - 1) / 2) { freq = (2 * PI * i) / (M * dx);       }
        else                  { freq = (2 * PI * (i - M)) / (M * dx); }
        // exponential of derivative operators
        exp_der[i] = cexp((-1) * Idt * a2 * freq * freq);
    }
    /* -------------------------------------------------------------------- */



    /* Apply Split step and solve separately nonlinear and linear part */
    /*******************************************************************/

    for (i = 0; i < N; i++)
    {
        // Apply exponential of potential together nonlinear part
        carrAbs2(M, S, abs2);
        rarrUpdate(M, V, inter, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);
        carrMultiply(M, stepexp, S, foward_fft);

        s = DftiComputeForward(desc, foward_fft);       // go to momentum space
        carrMultiply(M, exp_der, foward_fft, back_fft); // apply derivatives
        s = DftiComputeBackward(desc, back_fft);        // back to real space
        carrMultiply(M, stepexp, back_fft, Sstep);
        
        /* IMPROVEMENT
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the nonlinear part.
         *
         */

        // Trapezium rule in the nonlinear exponential part
        for (j = 0; j < M; j++)
        {
            abs2[j] += creal(Sstep[j]) * creal(Sstep[j]);
            abs2[j] += cimag(Sstep[j]) * cimag(Sstep[j]);
        }

        // factor / 4 due to trapezium rule
        rarrUpdate(M, V, inter / 4, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);
        carrMultiply(M, stepexp, S, foward_fft);

        s = DftiComputeForward(desc, foward_fft);       // go to momentum space
        carrMultiply(M, exp_der, foward_fft, back_fft); // apply derivatives
        s = DftiComputeBackward(desc, back_fft);        // back to real space
        carrMultiply(M, stepexp, back_fft, S);
    }

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(stepexp);
    free(foward_fft);
    free(back_fft);
    free(Sstep);
    free(abs2);
    free(out);
}




void GPCNLU(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname[], int n)
{
    FILE * out_data = fopen(fname, "w");
    
    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("ERROR: impossible to open file %s\n", fname);
        return;
    }

    unsigned int
        k,
        i,
        j;

    double complex
        Idt = 0.0 - dt * I;

    Carray Sstep = carrDef(M);

    // Used to apply nonlinear part
    Carray stepexp = carrDef(M);
    Carray linpart = carrDef(M);
    Rarray abs2    = rarrDef(M); // abs square of wave function

    // used to store matrix elements of linear part
    Carray upper = carrDef(M - 1);
    Carray lower = carrDef(M - 1);
    Carray mid   = carrDef(M - 1);
    Carray rhs   = carrDef(M - 1);



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
        carrAbs2(M, S, abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicLU(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

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

        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicLU(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        carrMultiply(M, linpart, stepexp, S);

        if (k == n) { record_step(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }
    }

    fclose(out_data);

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





void GPCNSM(int M, int N, double dx, double dt, double a2, double complex a1,
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
    Carray upper = carrDef(M - 1);
    Carray lower = carrDef(M - 1);
    Carray mid   = carrDef(M - 1);
    Carray rhs   = carrDef(M - 1);



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
        carrAbs2(M, S, abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

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

        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

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





void RK4DDT(int M, double inter, Carray S, Carray rhs)
{   // The Right-Hand-Side function to evolve with Runge-Kutta

    int i;

    for (i = 0; i < M; i++)
    {
        rhs[i] = inter * \
                 (creal(S[i]) * creal(S[i]) + cimag(S[i]) * cimag(S[i])) * S[i];
    }

    for (i = 0; i < M; i++) rhs[i] = (-1.0) * I * rhs[i];
}





void GPCNSMRK4_all(int M, int N, double dx, double dt, double a2,
     double complex a1, double inter, Rarray V, int cyclic, Cmatrix S)
{   // Evolve Gross-Pitaevskii using 4-th order Runge-Kutta
    // to deal with nonlinear part.
    int i,
        j;

    Carray karg = carrDef(M);
    Carray kans = carrDef(M);

    // Used to apply nonlinear part
    Carray linpart = carrDef(M);
    Rarray abs2    = rarrDef(M); // abs square of wave function

    // used to store matrix elements of linear part
    Carray upper = carrDef(M - 1);
    Carray lower = carrDef(M - 1);
    Carray mid   = carrDef(M - 1);
    Carray rhs   = carrDef(M - 1);



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



    dt = dt / 2; // Just needed in semi-step RK4 iteration
    for (i = 0; i < N; i++)
    {   // Apply exponential with nonlinear part

        // Compute k1 in kans
        RK4DDT(M, inter, S[i], kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k1
            linpart[j] = kans[j];
            // Argument to give k2
            karg[j] = S[i][j] + 0.5 * dt * kans[j];
        }
        
        // Compute k2 in kans
        RK4DDT(M, inter, karg, kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k2
            linpart[j] += 2 * kans[j];
            // Argument to give k3
            karg[j] = S[i][j] + 0.5 * dt * kans[j];
        }
        
        // Compute k3 in kans
        RK4DDT(M, inter, karg, kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k3
            linpart[j] += 2 * kans[j];
            // Argument to give k4
            karg[j] = S[i][j] + dt * kans[j];
        }

        // Compute k4 in kans
        RK4DDT(M, inter, karg, kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k4
            linpart[j] += kans[j];
        }
        
        for (j = 0; j < M; j++)
        {   // Finish RK4
            linpart[j] = S[i][j] + linpart[j] * dt / 6;
        }

        // Solve linear part
        CCSvec(M - 1, rhs_mat->vec, rhs_mat->col, rhs_mat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        // Compute k1 in kans
        RK4DDT(M, inter, linpart, kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k1
            S[i + 1][j] = kans[j];
            // Argument to give k2
            karg[j] = linpart[j] + 0.5 * dt * kans[j];
        }
        
        // Compute k2 in kans
        RK4DDT(M, inter, karg, kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k2
            S[i + 1][j] += 2 * kans[j];
            // Argument to give k3
            karg[j] = linpart[j] + 0.5 * dt * kans[j];
        }
        
        // Compute k3 in kans
        RK4DDT(M, inter, karg, kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k3
            S[i + 1][j] += 2 * kans[j];
            // Argument to give k4
            karg[j] = linpart[j] + dt * kans[j];
        }
        
        // Compute k4 in kans
        RK4DDT(M, inter, karg, kans);

        for (j = 0; j < M; j++)
        {   // Add up contribution from k2
            S[i + 1][j] += kans[j];
        }

        for (j = 0; j < M; j++)
        {   // Finish RK4
            S[i + 1][j] = linpart[j] + S[i + 1][j] * dt / 6;
        }
    }

    free(linpart);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    CCSFree(rhs_mat);

    free(karg);
    free(kans);
}
