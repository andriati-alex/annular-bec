#include "../include/GP_imagtime_integrator.h"





/* =======================================================================
 *
 *
 * MODULE OF FUNCTIONS THAT PROPAGATE THE GROSS-PITAEVSKII EQUATION
 * IN IMAGINARY TIME : t = - i T with T being real
 *
 *
 * i dS/dt = ( a2 (d^2/dx^2) + a1 (d/dx) + V(x) + inter |S|^2 ) S(x,t)
 *
 *
 * This module implement the mostly well know and used  integrators  to
 * the Gross-Pitaevskii(GP) equation. The routines are named as follows
 *
 *  - First 3 letters : IGP (from Gross-Pitaevskii)
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
 *  E end up with the energy on every time-step
 *
 *  CN methods supports both cyclic and zero boundary condition as
 *  identified by the cyclic(boolean) parameter.
 *
 * ======================================================================= */





void IGPFFT(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, Carray S, Carray E)
{   // Evolve the wave-function given an initial condition in S
    // on pure imaginary time to converge to an energy minimum.
    // Use FFT to compute derivatives on  linear  part  of  PDE
    // hence the boundary is required to be periodic



    int
        i,
        j,
        m = M - 1;
    
    
    
    MKL_LONG
        s;



    double
        freq,       // frequencies in Fourier space
        norm,       // initial norm
        NormStep,   // to renormalize on each time-step
        Idt = - dT; // factor to multiply on exponential after split-step



    double complex
        dt = - I  * dT; // pure imaginary time-step




    Rarray
        abs2 = rarrDef(M), // abs square of wave function
        out  = rarrDef(M); // hold linear and nonlinear potential



    Carray
        exp_der = carrDef(m),     // exponential of derivative operator
        stepexp = carrDef(M),     // Exponential of potential
        back_fft = carrDef(m),    // go back to position space
        forward_fft = carrDef(m); // go to frequency space



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
        exp_der[i] = cexp(Idt * a1 * freq * I - Idt * a2 * freq * freq);
    }
    /* ------------------------------------------------------------------- */



    /*   Apply Split step and solve separately nonlinear and linear part   */
    /*   ===============================================================   */

    for (i = 0; i < N; i++) {
        // Apply exponential ofone-body potential and nonlinear part
        rarrUpdate(M, V, inter, abs2, out);
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
    free(forward_fft);
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


    Carray
        // Exponential of nonlinear part
        stepexp = carrDef(M),
        // Hold the solution of linear system
        linpart = carrDef(M),
        // (cyclic)tridiagonal system from Crank-Nicolson
        upper = carrDef(M - 1),
        lower = carrDef(M - 1),
        mid   = carrDef(M - 1),
        // Vector of RHS of linear system
        rhs   = carrDef(M - 1);


    Rarray
        abs2 = rarrDef(M);


    CCSmat
        cnmat;



    /* Setup initial norm to be conserved and compute energy
     * ----------------------------------------------------- */
    carrAbs2(M, S, abs2);
    norm = sqrt(Rsimps(M, abs2, dx));
    E[0] = Functional(M, dx, a2, a1, inter / 2, V, S);
    /* ----------------------------------------------------- */



    // Configure the linear system from Crank-Nicolson scheme
    cnmat = CNmat(M, dx, dt, a2, a1, inter, V, cyclic, upper, lower, mid);



    for (i = 0; i < N; i++) {
        // Apply exponential with nonlinear part
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
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
    CCSFree(cnmat);
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

    Carray
        // Exponential of nonlinear part
        stepexp = carrDef(M),
        // Hold the solution of linear system
        linpart = carrDef(M),
        // (cyclic)tridiagonal system from Crank-Nicolson
        upper = carrDef(M - 1),
        lower = carrDef(M - 1),
        mid   = carrDef(M - 1),
        // Vector of RHS of linear system
        rhs   = carrDef(M - 1);


    Rarray
        abs2 = rarrDef(M);


    CCSmat
        cnmat;



    /* Setup initial norm to be conserved and compute energy
     * ----------------------------------------------------- */
    carrAbs2(M, S, abs2);
    norm = sqrt(Rsimps(M, abs2, dx));
    E[0] = Functional(M, dx, a2, a1, inter / 2, V, S);
    /* ----------------------------------------------------- */



    // Configure the linear system from Crank-Nicolson scheme
    cnmat = CNmat(M, dx, dt, a2, a1, inter, V, cyclic, upper, lower, mid);



    for (i = 0; i < N; i++) {
        // Apply exponential with nonlinear part
        rcarrExp(M, inter * Idt / 2, abs2, stepexp);
        carrMultiply(M, stepexp, S, linpart);

        // Solve linear part
        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
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
    CCSFree(cnmat);
}










/*  =====================================================================  *

                      NONLINEAR PART SOLVED BY RUNGE-KUTTA

    =====================================================================  */

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



void NonLinearVIDDT(int M, double t, Carray Psi, Carray FullPot, Carray Dpsi)
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
        Dpsi[i] = - (FullPot[0] * mod_squared + FullPot[i + 1] ) * Psi[i];
    }
}










void IGPCNSMRK4(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, int cyclic, Carray S, Carray E)
{   // Evolve Gross-Pitaevskii using 4-th order Runge-Kutta
    // to deal with nonlinear  part.  Solve  Crank-Nicolson
    // linear system with Sherman-Morrison formula


    int i,
        j;


    double
        norm,     // initial norm to be kept constant
        NormStep; // renorm each time step


    double complex
        dt,
        interv[1];


    Rarray
        abs2 = rarrDef(M); // abs square of wave function
    

    Carray
        // hold solution of linear part
        linpart = carrDef(M),
        // (cyclic)tridiagonal matrix from Crank-Nicolson
        upper   = carrDef(M - 1),
        lower   = carrDef(M - 1),
        mid     = carrDef(M - 1),
        // RHS of linear system at each time step
        rhs     = carrDef(M - 1);

    CCSmat
        cnmat;


    dt = - I * dT;     // Set time to pure imaginary
    interv[0] = inter; // setup extra arguments to RK4 derivative


    
    /* Setup initial norm to be conserved and compute energy
     * ----------------------------------------------------- */
    carrAbs2(M, S, abs2);
    norm = sqrt(Rsimps(M, abs2, dx));
    E[0] = Functional(M, dx, a2, a1, inter / 2, V, S);
    /* ----------------------------------------------------- */



    // Configure the linear system from Crank-Nicolson scheme
    cnmat = CNmat(M, dx, dt, a2, a1, inter, V, cyclic, upper, lower, mid);



    for (i = 0; i < N; i++)
    {
        RK4step(M, dT/2, 0, S, interv, linpart, NonLinearIDDT);
        
        // Solve linear part (nabla ^ 2 part)
        CCSvec(M - 1, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
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

        if ((i + 1) % 2000 == 0 && i < 10000)
        {
            dt = dt * (1 + 0.2);
            dT = dT * (1 + 0.2);
            CCSFree(cnmat); // Erase old matrix to setup new one
            cnmat = CNmat( M, dx, dt, a2, a1, inter,
                    V, cyclic, upper, lower, mid);
        }
    }

    free(linpart);
    free(upper);
    free(lower);
    free(abs2);
    free(mid);
    free(rhs);
    CCSFree(cnmat);
}










void IGPFFTRK4(int M, int N, double dx, double dT, double a2, double complex a1,
     double inter, Rarray V, Carray S, Carray E)
{   // Evolve the wave-function given an initial condition in S
    // on pure imaginary time to converge to an energy minimum.
    // Use FFT to compute derivatives on  linear  part  of  PDE
    // hence the boundary is required to be periodic. Nonlinear
    // part is solved by 4th order Runge-Kutta



    int
        i,
        j,
        m = M - 1;
    
    
    
    MKL_LONG
        s;



    double
        freq,       // frequencies in Fourier space
        norm,       // initial norm
        NormStep,   // to renormalize on each time-step
        Idt = - dT; // factor to multiply on exponential after split-step



    double complex
        dt = - I * dT;



    Rarray
        abs2 = rarrDef(M); // abs square of wave function



    Carray
        argRK4  = carrDef(M),     // hold linear and nonlinear potential
        exp_der = carrDef(m),     // exponential of derivative operator
        back_fft = carrDef(m),    // go back to position space
        forward_fft = carrDef(m), // go to frequency space
        FullPot = carrDef(M + 1);



    // Setup RHS of potential splitted-step part of derivatives
    FullPot[0] = inter;
    for (i = 0; i < M; i++) FullPot[i + 1] = V[i];



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
        exp_der[i] = cexp(Idt * a1 * freq * I - Idt * a2 * freq * freq);
    }
    /* ------------------------------------------------------------------- */



    /*   Apply Split step and solve separately nonlinear and linear part   */
    /*   ===============================================================   */

    for (i = 0; i < N; i++)
    {
        RK4step(M, dT/2, 0, S, FullPot, argRK4, NonLinearVIDDT);
        carrCopy(m, argRK4, forward_fft);

        // go to momentum space
        s = DftiComputeForward(desc, forward_fft);
        // apply exponential of derivatives
        carrMultiply(m, exp_der, forward_fft, back_fft);
        // go back to position space
        s = DftiComputeBackward(desc, back_fft);
        carrCopy(m, back_fft, argRK4);
        argRK4[m] = argRK4[0];
        
        RK4step(M, dT/2, 0, argRK4, FullPot, S, NonLinearVIDDT);

        carrAbs2(M, S, abs2);

        // Renormalization
        NormStep = norm / sqrt(Rsimps(M, abs2, dx));
        for (j = 0; j < M; j++) S[j] = NormStep * S[j];

        // Energy
        E[i + 1] = Functional(M, dx, a2, a1, inter / 2, V, S);
    }

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(forward_fft);
    free(back_fft);
    free(abs2);
    free(argRK4);
    free(FullPot);
}
