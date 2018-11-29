#include "../include/GP_realtime_integrator.h"





/* =======================================================================
 *
 *
 * MODULE OF FUNCTIONS THAT INTEGRATE THE GROSS-PITAEVSKII EQUATION
 *
 *
 * i dS/dt = ( a2 (d^2/dx^2) + a1 (d/dx) + V(x) + inter |S|^2 ) S(x,t)
 *
 *
 * This module implement the mostly well know and used  integrators  to
 * the Gross-Pitaevskii(GP) equation. The routines are named as follows
 *
 *  - First 2 letters : GP (from Gross-Pitaevskii)
 *
 *  - Next letters : identify the methods
 *      (1) CN for Crank-Nicolson finite differences scheme
 *      (2) FFT for Fast Fourier transform to deal with derivatives
 *
 *  - The next letters apply to CN on how to solve the linear system
 *      (1) SM for Sherman-Morrison formula.
 *      (2) LU for the decomposition.
 *  
 *  - For those who use 4th order Runge-Kutta method to integrate the
 *    nonlinear part from the split-step, we have additionally RK4 as
 *    a suffix, otherwise,  it  is  done going one step forward using
 *    integration by rectangles and doing again the same steps  using
 *    integration by trapezium rule
 *
 *  M is the number of discretized points (size of arrays)
 *  N is the number of time-steps to be propagated the initial condition
 *
 *  CN methods supports both cyclic and zero boundary condition as
 *  identified by the cyclic(boolean) parameter.
 *
 * ======================================================================= */




void record_step(FILE * f, int M, Carray v)
{

//  Given an open file f write the complex array v in a line

    int j;

    for (j = 0; j < M; j ++)
    {   
        // write in a format suitable to import with numpy
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





void GPFFT(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, Carray S, char fname[], int n)
{

//  Evolve the wave-function given an initial condition in S
//  that is overwritten at each time-step.  The  results are
//  recorded in a file named 'fname' on every 'n' steps. Use
//  FFT to compute linear derivatives part of  PDE hence the
//  boundary is consider to be periodic in the last position


    // File to write every n step the time-step solution
    FILE * out_data = fopen(fname, "w");

    if (out_data == NULL)
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    record_step(out_data, M, S);


    // Reader of screen printing
    printf("\n\n\n");
    printf("     step            Energy                   Norm");
    SepLine();


    int
        k,
        i,
        j,
        m = M - 1;


    MKL_LONG
        s; // status of called MKL FFT functions


    double
        freq;


    double complex
        E,
        Idt = 0.0 - dt * I;



    Rarray
        abs2 = rarrDef(M), // abs square of wave function
        out  = rarrDef(M); // output of linear + nonlinear potential



    Carray
        exp_der = carrDef(m),     // exponential of derivative operator
        stepexp = carrDef(M),     // Exponential of potential
        forward_fft = carrDef(m), // go to frequency space
        back_fft = carrDef(m),    // back to position space
        Sstep = carrDef(M);       // hold one step to integrate by trapezium



    /* setup descriptor (MKL implementation of FFT)
     * -------------------------------------------------------------------- */
    DFTI_DESCRIPTOR_HANDLE desc;
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, m);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(m));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(m));
    s = DftiCommitDescriptor(desc);
    /* -------------------------------------------------------------------- */



    /* setup Fourier Frequencies and the exponential of derivative operator
     * -------------------------------------------------------------------- */
    for (i = 0; i < m; i++) {
        if (i <= (m - 1) / 2) { freq = (2 * PI * i) / (m * dx);       }
        else                  { freq = (2 * PI * (i - m)) / (m * dx); }
        // exponential of derivative operators
        exp_der[i] = cexp(Idt * a1 * freq * I - Idt * a2 * freq * freq);
    }
    /* -------------------------------------------------------------------- */



    /*  Apply Split step and solve separately nonlinear and linear part  *
     *  ===============================================================  */



    k = 1;
    for (i = 0; i < N; i++)
    {
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Functional(M, dx, a2, a1, inter / 2, V, S);
        printf(" \n  %7d          ", i);
        printf("%15.7E          ", creal(E));
        printf("%15.7E          ", Rsimps(M, abs2, dx));



        // Apply exponential of potential together nonlinear part
        rarrUpdate(M, V, inter, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);
        carrMultiply(m, stepexp, S, forward_fft);

        // go to momentum space
        s = DftiComputeForward(desc, forward_fft);
        // apply exponential of derivatives
        carrMultiply(m, exp_der, forward_fft, back_fft);
        // go back to position space
        s = DftiComputeBackward(desc, back_fft);

        carrCopy(m, back_fft, Sstep);
        Sstep[m] = Sstep[0];

        carrAbs2(M, Sstep, abs2);
        rarrUpdate(M, V, inter, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);

        carrMultiply(m, stepexp, back_fft, Sstep);
        Sstep[m] = Sstep[0];



        /* IMPROVEMENT
         * -----------
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the nonlinear  part.
         *
         * ------------------------------------------------------------ */



        // Trapezium rule in the nonlinear exponential part
        carrAbs2(M, S, abs2);
        for (j = 0; j < M; j++)
        {
            abs2[j] += creal(Sstep[j]) * creal(Sstep[j]);
            abs2[j] += cimag(Sstep[j]) * cimag(Sstep[j]);
        }

        // factor / 2 due to trapezium rule
        rarrUpdate(M, V, inter / 2, abs2, out);
        rcarrExp(M, Idt / 2, out, stepexp);
        carrMultiply(m, stepexp, S, forward_fft);

        // go to momentum space
        s = DftiComputeForward(desc, forward_fft);
        // apply exponential of derivatives
        carrMultiply(m, exp_der, forward_fft, back_fft);
        // go back to position space
        s = DftiComputeBackward(desc, back_fft);
        carrMultiply(m, stepexp, back_fft, S);
        S[m] = S[0];



        // RECORD solution

        if (k == n) { record_step(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }
    }

    carrAbs2(M, S, abs2);
    E = Functional(M, dx, a2, a1, inter / 2, V, S);
    printf(" \n  %7d          ", N);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));
    
    SepLine();

    fclose(out_data);

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(stepexp);
    free(forward_fft);
    free(back_fft);
    free(Sstep);
    free(abs2);
    free(out);
}





void GPCNLU(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname[], int n)
{

//  Evolve the wave-function given an initial condition  in S
//  that is overwritten at each time-step.  The  results  are
//  recorded in a file named 'fname' on every 'n' steps.  Use
//  Crank-Nicolson discretization scheme to  compute   linear
//  part of potential and derivatives of the PDE. 'Cyclic' is
//  a boolean argument define  whether  the  boundary is zero
//  or periodic on last position point.


    // File to write every n step the time-step solution
    FILE * out_data = fopen(fname, "w");



    if (out_data == NULL)
    {
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    record_step(out_data, M, S);
    
    
    
    // Reader of screen printing
    printf("\n\n\n");
    printf("     step            Energy                   Norm");
    SepLine();


    unsigned int
        k,
        i,
        j;



    // Factor to multiply in nonlinear exponential part after split-step
    double complex
        E,
        Idt = 0.0 - dt * I;



    Rarray
        abs2 = rarrDef(M); // abs square of wave function



    Carray
        // go back and integrate nonlinear part by trapezium
        Sstep = carrDef(M),
        // Used to apply nonlinear part
        stepexp = carrDef(M),
        // hold solution of linear system
        linpart = carrDef(M),
        // (cyclic)tridiagonal matrix from Crank-Nicolson
        upper = carrDef(M - 1),
        lower = carrDef(M - 1),
        mid   = carrDef(M - 1),
        // RHS of linear system to solve
        rhs   = carrDef(M - 1);
    
    
    
    CCSmat
        cnmat;



            /*  ==========================================

                       Setup Right-Hand-Side matrix

                ==========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M - 1, upper, dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                        }

    // fill lower diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                        }

    // Store in CCS format
    cnmat = CyclicToCCS(M - 1, upper, lower, mid);



            /*  ===========================================

                      Setup Cyclic tridiagonal matrix

                ===========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(M - 1, upper, -dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                          }

    // fill lower diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                          }



    /*  Apply Split step and solve separately nonlinear and linear part  *
     *  ===============================================================  */



    k = 1;
    for (i = 0; i < N; i++)
    {
        
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Functional(M, dx, a2, a1, inter / 2, V, S);
        printf(" \n  %7d          ", i);
        printf("%15.7E          ", creal(E));
        printf("%15.7E          ", Rsimps(M, abs2, dx));



        // Apply exponential with nonlinear part
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicLU(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary
        
        // Apply exponential with nonlinear part again
        carrAbs2(M, linpart, abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, linpart, Sstep);



        /* IMPROVEMENT
         * -----------
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the  nonlinear part
         *
         * ------------------------------------------------------------ */



        // Trapezium rule in the nonlinear exponential part
        carrAbs2(M, S, abs2);
        for (j = 0; j < M; j++)
        {
            abs2[j] += creal(Sstep[j]) * creal(Sstep[j]);
            abs2[j] += cimag(Sstep[j]) * cimag(Sstep[j]);
        }

        // Extra factor 1/2 due to trapezium rule
        rcarrExp(M, inter * Idt / 4, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicLU(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        carrMultiply(M, linpart, stepexp, S);

        // record data every n steps
        if (k == n) { record_step(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }

    }

    carrAbs2(M, S, abs2);
    E = Functional(M, dx, a2, a1, inter / 2, V, S);
    printf(" \n  %7d          ", N);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));
    
    SepLine();

    fclose(out_data);

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    free(Sstep);
    CCSFree(cnmat);
}





void GPCNSM(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname [], int n)
{

//  Evolve the wave-function given an initial condition in S
//  that is overwritten at each time-step.  The  results are
//  recorded in a file named 'fname' on every 'n' steps. Use
//  Crank-Nicolson discretization scheme to  compute  linear
//  part of potential and derivatives of the PDE.   'Cyclic'
//  is a boolean argument define  whether  the  boundary  is
//  zero or periodic on last position point.


    // File to write every n step the time-step solution
    FILE * out_data = fopen(fname, "w");



    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    record_step(out_data, M, S);
    
    
    
    // Reader of screen printing
    printf("\n\n\n");
    printf("     step            Energy                   Norm");
    SepLine();



    unsigned int
        k,
        i,
        j;



    // Factor to multiply in nonlinear exponential part after split-step
    double complex
        E,
        Idt = 0.0 - dt * I;



    Rarray
        abs2 = rarrDef(M); // abs square of wave function



    Carray
        // go back and integrate nonlinear part by trapezium
        Sstep = carrDef(M),
        // Used to apply nonlinear part
        stepexp = carrDef(M),
        // hold solution of linear system
        linpart = carrDef(M),
        // (cyclic)tridiagonal matrix from Crank-Nicolson
        upper = carrDef(M - 1),
        lower = carrDef(M - 1),
        mid   = carrDef(M - 1),
        // RHS of linear system to solve
        rhs   = carrDef(M - 1);
    
    
    
    CCSmat
        cnmat;



            /*  ==========================================

                       Setup Right-Hand-Side matrix

                ==========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M - 1, upper, dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                        }

    // fill lower diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                        }

    // Store in CCS format
    cnmat = CyclicToCCS(M - 1, upper, lower, mid);



            /*  ===========================================

                      Setup Cyclic tridiagonal matrix

                ===========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(M - 1, upper, -dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                          }

    // fill lower diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                          }


    
    /*  Apply Split step and solve separately nonlinear and linear part  *
     *  ===============================================================  */



    k = 1;
    for (i = 0; i < N; i++)
    {

        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Functional(M, dx, a2, a1, inter / 2, V, S);
        printf(" \n  %7d          ", i);
        printf("%15.7E          ", creal(E));
        printf("%15.7E          ", Rsimps(M, abs2, dx));



        // Apply exponential with nonlinear part
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        // Apply exponential with nonlinear part again
        carrAbs2(M, linpart, abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, linpart, stepexp, Sstep);



        /* IMPROVEMENT
         * -----------
         *
         * The steps above use rectangular integration of the nonlinear
         * part,  thus with an approximate solution to the next step we
         * can improve by using trapezium method in the nonlinear part.
         *
         * ------------------------------------------------------------ */



        carrAbs2(M, S, abs2);
        // Trapezium rule in the nonlinear exponential part
        for (j = 0; j < M; j++)
        {
            abs2[j] += creal(Sstep[j]) * creal(Sstep[j]);
            abs2[j] += cimag(Sstep[j]) * cimag(Sstep[j]);
        }

        // Extra factor 1/2 due to trapezium rule
        rcarrExp(M, inter * Idt / 4, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        carrMultiply(M, linpart, stepexp, S);

        // record data every n steps
        if (k == n) { record_step(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }
    }

    carrAbs2(M, S, abs2);
    E = Functional(M, dx, a2, a1, inter / 2, V, S);
    printf(" \n  %7d          ", N);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));
    
    SepLine();

    fclose(out_data);

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    free(Sstep);
    CCSFree(cnmat);
}





/*  =====================================================================  *

                      NONLINEAR PART SOLVED BY RUNGE-KUTTA

    =====================================================================  */





void NonLinearDDT(int M, double t, Carray Psi, Carray inter, Carray Dpsi)
{   // The Right-Hand-Side of derivative of non-linear part is computed
    // on Dpsi after used split-step technique. A function  to  be used
    // on 4th order Runge-Kutta routine. The array of extra  parameters
    // contains only one inter[0] with contact interaction strength.

    int i;

    // sum real and imaginary part squared to get the modulus squared
    double mod_squared;

    for (i = 0; i < M; i++)
    {
        mod_squared  = creal(Psi[i]) * creal(Psi[i]);
        mod_squared += cimag(Psi[i]) * cimag(Psi[i]);
        Dpsi[i] = - I * inter[0] * mod_squared * Psi[i];
    }
}





void NonLinearVDDT(int M, double t, Carray Psi, Carray FullPot, Carray Dpsi)
{   // The Right-Hand-Side of derivative of non-linear and linear  potential
    // part. As a extra argument take the interaction strength in FullPot[0]
    // and the linear potential computed at position i in FullPot[i + 1]

    int i;

    // sum real and imaginary part squared to get the modulus squared
    double mod_squared;

    for (i = 0; i < M; i++)
    {
        mod_squared  = creal(Psi[i]) * creal(Psi[i]);
        mod_squared += cimag(Psi[i]) * cimag(Psi[i]);
        Dpsi[i] = - I * (FullPot[0] * mod_squared + FullPot[i + 1] ) * Psi[i];
    }
}





void GPCNSMRK4(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname [], int n)
{   // Evolve Gross-Pitaevskii using 4-th order Runge-Kutta to deal
    // with nonlinear part and Crank-Nicolson  to  the  linear part


    // File to write every n step the time-step solution
    FILE * out_data = fopen(fname, "w");



    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    record_step(out_data, M, S);



    // Reader of screen printing
    printf("\n\n\n");
    printf("     step            Energy                   Norm");
    SepLine();



    int k,
        i,
        j;



    double complex
        E,
        interv[1];



    Rarray
        abs2 = rarrDef(M);



    Carray
        linpart = carrDef(M),
        // (cyclic)tridiagonal system
        upper = carrDef(M - 1),
        lower = carrDef(M - 1),
        mid   = carrDef(M - 1),
        // RHS to slve the linear system at each time step
        rhs   = carrDef(M - 1);



    CCSmat
        cnmat;



    interv[0] = inter; // extra arg to RK4 derivatives



            /*  ==========================================

                       Setup Right-Hand-Side matrix

                ==========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M - 1, upper, dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                        }

    // fill lower diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                        }

    // Store in CCS format
    cnmat = CyclicToCCS(M - 1, upper, lower, mid);



            /*  ===========================================

                      Setup Cyclic tridiagonal matrix

                ===========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(M - 1, upper, -dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                          }

    // fill lower diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                          }



    /*  Apply Split step and solve separately nonlinear and linear part  *
     *  ===============================================================  */



    k = 1;
    for (i = 0; i < N; i++)
    {
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Functional(M, dx, a2, a1, inter / 2, V, S);
        printf(" \n  %7d          ", i);
        printf("%15.7E          ", creal(E));
        printf("%15.7E          ", Rsimps(M, abs2, dx));



        RK4step(M, dt/2, 0, S, interv, linpart, NonLinearDDT);
        
        // Solve linear part (nabla ^ 2 part + onebody potential)
        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        RK4step(M, dt/2, 0, linpart, interv, S, NonLinearDDT);



        // record data every n steps
        if (k == n) { record_step(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }
    }

    carrAbs2(M, S, abs2);
    E = Functional(M, dx, a2, a1, inter / 2, V, S);
    printf(" \n  %7d          ", N);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));

    SepLine();

    fclose(out_data);

    free(linpart);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    free(abs2);
    CCSFree(cnmat);
}





void GPFFTRK4(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, Carray S, char fname[], int n)
{   // Evolve the wave-function given an initial condition in S
    // that is overwritten at each time-step.  The  results are
    // recorded in a file named 'fname' on every 'n' steps. Use
    // FFT to compute linear derivatives part of  PDE hence the
    // boundary is consider to be periodic in the last position
    // Solve nonlinear part using 4th order Runge-Kutta



    // File to write every n step the time-step solution
    FILE * out_data = fopen(fname, "w");



    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    record_step(out_data, M, S);



    // Reader of screen printing
    printf("\n\n\n");
    printf("     step            Energy                   Norm");
    SepLine();



    int
        k,
        i,
        j,
        m = M - 1;



    MKL_LONG
        s;



    double
        freq;



    double complex
        E,
        Idt = 0 - dt * I;


    Rarray
        abs2 = rarrDef(M);



    Carray
        argRK4  = carrDef(M),     // output of linear and nonlinear pot
        exp_der = carrDef(m),     // exponential of derivative operator
        forward_fft = carrDef(m), // go to frequency space
        back_fft = carrDef(m),    // back to position space
        FullPot = carrDef(M + 1);



    // Setup RHS of potential for nonlinear + trap potential
    FullPot[0] = inter;
    for (i = 0; i < M; i++) FullPot[i + 1] = V[i];
    


    /* setup descriptor (MKL implementation of FFT)
     * -------------------------------------------------------------------- */
    DFTI_DESCRIPTOR_HANDLE desc;
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, m);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(m));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(m));
    s = DftiCommitDescriptor(desc);
    /* -------------------------------------------------------------------- */



    /* setup Fourier Frequencies and the exponential of derivative operator
     * -------------------------------------------------------------------- */
    for (i = 0; i < m; i++) {
        if (i <= (m - 1) / 2) { freq = (2 * PI * i) / (m * dx);       }
        else                  { freq = (2 * PI * (i - m)) / (m * dx); }
        // exponential of derivative operators
        exp_der[i] = cexp(Idt * a1 * freq * I - Idt * a2 * freq * freq);
    }
    /* -------------------------------------------------------------------- */



    /*  Apply Split step and solve separately nonlinear and linear part  *
     *  ===============================================================  */



    k = 1;
    for (i = 0; i < N; i++)
    {

        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Functional(M, dx, a2, a1, inter / 2, V, S);
        printf(" \n  %7d          ", i);
        printf("%15.7E          ", creal(E));
        printf("%15.7E          ", Rsimps(M, abs2, dx));



        RK4step(M, dt/2, 0, S, FullPot, argRK4, NonLinearVDDT);
        carrCopy(m, argRK4, forward_fft);

        // go to momentum space
        s = DftiComputeForward(desc, forward_fft);
        // apply exponential of derivatives
        carrMultiply(m, exp_der, forward_fft, back_fft);
        // go back to position space
        s = DftiComputeBackward(desc, back_fft);
        carrCopy(m, back_fft, argRK4);
        argRK4[m] = argRK4[0]; // cyclic condition

        RK4step(M, dt/2, 0, argRK4, FullPot, S, NonLinearVDDT);

        // record data every n steps
        if (k == n) { record_step(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }

    }

    carrAbs2(M, S, abs2);
    E = Functional(M, dx, a2, a1, inter / 2, V, S);
    printf(" \n  %7d          ", N);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));

    SepLine();

    fclose(out_data);

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(forward_fft);
    free(back_fft);
    free(argRK4);
    free(abs2);
    free(FullPot);
}










void CFDS(int M, int N, double dx, double dt, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, char fname [], int n)
{

//  Evolve the wave-function given an initial condition in S
//  that is overwritten at each time-step.  The  results are
//  recorded in a file named 'fname' on every 'n' steps. Use
//  Crank-Nicolson discretization scheme to  compute  linear
//  part of potential and derivatives of the PDE.   'Cyclic'
//  is a boolean argument define  whether  the  boundary  is
//  zero or periodic on last position point.


    // File to write every n step the time-step solution
    FILE * out_data = fopen(fname, "w");



    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    record_step(out_data, M, S);



    // Reader of screen printing
    printf("\n\n\n");
    printf("     step            Energy                   Norm");
    SepLine();



    unsigned int
        k,
        i,
        j,
        iter,
        condition;


    double
        tol;



    // Factor to multiply in nonlinear exponential part after split-step
    double complex
        aux,
        Idt = 0.0 - dt * I;



    Rarray
        abs2 = rarrDef(M); // abs square of wave function



    Carray
        // go back and integrate nonlinear part by trapezium
        Sstep = carrDef(M),
        // Used to apply nonlinear part
        stepexp = carrDef(M),
        // hold solution of linear system
        linpart = carrDef(M),
        // (cyclic)tridiagonal matrix from Crank-Nicolson
        upper = carrDef(M - 1),
        lower = carrDef(M - 1),
        mid   = carrDef(M - 1),
        // RHS of linear system to solve
        rhs   = carrDef(M - 1);



    CCSmat
        cnmat;



            /*  ==========================================

                       Setup Right-Hand-Side matrix

                ==========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(M - 1, upper, dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                        }

    // fill lower diagonal
    carrFill(M - 1, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                        }

    // Store in CCS format
    cnmat = CyclicToCCS(M - 1, upper, lower, mid);



            /*  ===========================================

                      Setup Cyclic tridiagonal matrix

                ===========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(M - 1, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(M - 1, upper, -dt / 2, V, mid);

    // fill upper diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[M-2] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[M-2] = 0;                                          }

    // fill lower diagonal
    carrFill(M - 1, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[M-2] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[M-2] = 0;                                          }


    
    /*  Apply Split step and solve separately nonlinear and linear part  *
     *  ===============================================================  */



    k = 1;
    for (i = 0; i < N; i++)
    {

        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        aux = Functional(M, dx, a2, a1, inter / 2, V, S);
        printf(" \n  %7d          ", i);
        printf("%15.7E          ", creal(aux));
        printf("%15.7E          ", Rsimps(M, abs2, dx));



        // Apply exponential with nonlinear part
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicSM(M - 1, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        // Apply exponential with nonlinear part again
        carrAbs2(M, linpart, abs2);
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, linpart, stepexp, Sstep);



        // Change main diagonal to do iterations
        for (j = 0; j < M - 1; j ++)
        {
            mid[j] = mid[j] - 0.25 * dt * inter * cabs(S[j]) * cabs(S[j]);
        }


        condition = 1;
        iter = 0;
        while (condition)
        {

            carrCopy(M, Sstep, linpart);
            CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, S, rhs);

            for (j = 0; j  < M - 1; j++)
            {
                aux = cabs(Sstep[j]) * cabs(Sstep[j]) * (Sstep[j] + S[j]);
                aux = aux + cabs(S[j]) * cabs(S[j]) * S[j];
                rhs[j] = rhs[j] + 0.25 * dt * inter * aux;
            }

            triCyclicSM(M - 1, upper, lower, mid, rhs, Sstep);
            Sstep[M-1] = Sstep[0];

            condition = 0;
            for (j = 0; j  < M; j++)
            {
                tol = 1E-9 * cabs(linpart[j]) + 1E-12;
                if (cabs(Sstep[j] - linpart[j]) > tol) condition = 1;
            }

            iter = iter + 1;

        }

        // correct the main diagonal
        for (j = 0; j < M - 1; j ++)
        {
            mid[j] = I + a2 * dt / dx / dx - dt * V[j] / 2;
        }

        carrCopy(M, Sstep, S);
        
        // record data every n steps
        if (k == n) { record_step(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }

    }

    carrAbs2(M, S, abs2);
    aux = Functional(M, dx, a2, a1, inter / 2, V, S);
    printf(" \n  %7d          ", N);
    printf("%15.7E          ", creal(aux));
    printf("%15.7E          ", Rsimps(M, abs2, dx));

    SepLine();

    fclose(out_data);

    free(stepexp);
    free(linpart);
    free(abs2);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    free(Sstep);
    CCSFree(cnmat);
}
