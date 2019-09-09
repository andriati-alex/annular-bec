#include "realtime_integrator.h"





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





void SSFFT(EqDataPkg EQ, int N, double dt, Carray S, char fname[], int n)
{

/** Evolve the wave-function given an initial condition in S that  is
  * overwritten at each time-step. The results are recorded in a file
  * named 'fname' on every 'n' steps.
  *
  * The numerical method is a Split-Step which separate the linear operator
  * involving derivatives from the potential part, this last one including
  * a linear and nonlinear part. It is used FFT to compute the derivarives **/



    int
        k,
        i,
        j,
        M,
        m;

    MKL_LONG
        s;

    double
        a2,
        dx,
        g,
        freq;

    double complex
        E,
        a1,
        Idt = 0.0 - dt * I;

    DFTI_DESCRIPTOR_HANDLE
        desc;

    Rarray
        V,
        abs2,
        step_pot;

    Carray
        exp_der,
        exp_pot,
        forward_fft,
        back_fft;

    FILE
        * out_data = fopen(fname, "w");



    M = EQ->Mpos;   // grid size including boudaries
    m = M - 1;      // grid size excluding boudaries

    abs2 = rarrDef(M);      // abs square of wave function
    step_pot  = rarrDef(M); // time-dependent potential (nonlinear+linear)

    exp_der = carrDef(m);     // Exponential of derivative operators
    exp_pot = carrDef(M);     // Exponential of potential
    forward_fft = carrDef(m);
    back_fft = carrDef(m);

    // Open file to write solution at every n time steps

    out_data = fopen(fname, "w");

    if (out_data == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line

    carr_inline(out_data, M, S);

    // unpack equation parameters from structure

    a2 = EQ->a2;
    a1 = EQ->a1;
    dx = EQ->dx;
    g = EQ->inter;
    V = EQ->V;



    // setup descriptor (MKL implementation of FFT)
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, m);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(m));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(m));
    s = DftiCommitDescriptor(desc);



    // setup Fourier Frequencies and the exponential of derivative operator
    // The convention for ordering the frequencies is starting from zero go
    // over positive frequencies and after reached the max. start from  the
    // max. negative frequency going back to zero.
    for (i = 0; i < m; i++)
    {
        if (i <= (m - 1) / 2) { freq = (2 * PI * i) / (m * dx);       }
        else                  { freq = (2 * PI * (i - m)) / (m * dx); }
        // exponential of derivative operators
        exp_der[i] = cexp(Idt * a1 * freq * I - Idt * a2 * freq * freq);
    }



    // Header of screen printing
    printf("\n\n\n");
    printf("     time            Energy                   Norm");
    sepline();



    k = 1;
    for (i = 0; i < N; i++)
    {
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control

        E = Energy(M, dx, a2, a1, g, V, S);
        if ( i % 50 == 0 )
        {
            printf(" \n  %.4lf          ", i*dt);
            printf("%15.7E          ", creal(E));
            printf("%15.7E          ", Rsimps(M, abs2, dx));
        }



        // Apply exponential of potential part (linear and nonlinear)
        // When copying data to use Fourier transform it is not used
        // the boundary grid point assumed to be periodic
        rarrUpdate(M, V, g, abs2, step_pot);
        rcarrExp(M, Idt / 2, step_pot, exp_pot);
        carrMultiply(m, exp_pot, S, forward_fft);



        // go to momentum space
        s = DftiComputeForward(desc, forward_fft);
        // apply exponential of derivatives
        carrMultiply(m, exp_der, forward_fft, back_fft);
        // go back to position space
        s = DftiComputeBackward(desc, back_fft);
        carrCopy(m, back_fft, S);
        S[m] = S[0]; //boundary



        // Apply again the full potential part
        carrAbs2(M, S, abs2);
        rarrUpdate(M, V, g, abs2, step_pot);
        rcarrExp(M, Idt / 2, step_pot, exp_pot);
        carrMultiply(m, exp_pot, back_fft, S);
        S[m] = S[0]; //boundary



        // RECORD solution if required
        if (k == n) { carr_inline(out_data, M, S); k = 1; }
        else        { k = k + 1; }
    }

    carrAbs2(M, S, abs2);
    E = Energy(M, dx, a2, a1, g, V, S);
    printf(" \n  %.4lf          ", N*dt);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));

    sepline();

    fclose(out_data);

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(exp_pot);
    free(forward_fft);
    free(back_fft);
    free(abs2);
    free(step_pot);
}





void SSCNLU(EqDataPkg EQ, int N, double dt, int cyclic, Carray S,
     char fname[], int n)
{

/** Evolve the wave-function given an initial condition in S that  is
  * overwritten at each time-step. The results are recorded in a file
  * named 'fname' on every 'n' steps.
  *
  * The numerical method is a Split-Step which separate the linear operator
  * from the nonlinear part. It is used Crank-Nicolson  finite  differences
  * scheme to transform all the linear part in a linear system of equations
  * The linear system of equations is solved using LU decompodition  for  a
  * tridiagonal system that may include cyclic terms from periodic boundary
  * conditions. For more information on this see
  *
  * Thiab R. Taha, "Solution of Periodic Linear Systems on a Hypercube",
  * IEEE, doi:10.1109/DMCC.1990.555404                               **/



    unsigned int
        k,
        i,
        M,
        m,
        j;

    double
        a2,
        dx,
        g;

    double complex
        E,
        a1,
        Idt;

    Rarray
        V,
        abs2;

    Carray
        exp_pot,
        linpart,
        upper,
        lower,
        mid,
        rhs;

    CCSmat
        cnmat;

    FILE
        * out_data;



    M = EQ->Mpos; // grid size
    m = M - 1;    // grid size excluding the boundaries

    abs2 = rarrDef(M);    // squared modulus of wave function at grid points
    exp_pot = carrDef(M); // nonlinear potential part exponential in SS
    linpart = carrDef(M); // linear part solution by Finite Differences

    // tridiagonal system with additional  cyclic  terms depending on
    // the boundary conditions given.  In  positive case for periodic
    // boundary use the last element of 'upper' and 'lower' diagonals
    // to store the last column first line element in 'upper' and the
    // last line first column element in 'lower'.
    upper = carrDef(m);
    lower = carrDef(m);
    mid   = carrDef(m);

    rhs   = carrDef(m); // RHS of linear system to solve

    out_data = fopen(fname, "w");
    if (out_data == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    carr_inline(out_data, M, S);

    // unpack equation parameters from structure
    a2 = EQ->a2;
    a1 = EQ->a1;
    dx = EQ->dx;
    g = EQ->inter;
    V = EQ->V;
    Idt = 0.0 - dt * I;



            /*  ===========================================

                Setup RHS matrix from Crank-Nicolson Scheme

                ===========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(m, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(m, upper, dt / 2, V, mid);

    // fill upper diagonal
    carrFill(m, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[m-1] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[m-1] = 0;                                        }

    // fill lower diagonal
    carrFill(m, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[m-1] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[m-1] = 0;                                        }

    // Store in Compressed-Comlumn-Storage format
    cnmat = cyclic2CCS(m, upper, lower, mid);



          /*  ================================================

              Setup Cyclic tridiagonal matrix of linear system

              ================================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(m, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(m, upper, -dt / 2, V, mid);

    // fill upper diagonal
    carrFill(m, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[m-1] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[m-1] = 0;                                          }

    // fill lower diagonal
    carrFill(m, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[m-1] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[m-1] = 0;                                          }



    // Header of the screen output
    printf("\n\n\n");
    printf("     Time            Energy                   Norm");
    sepline();



    k = 1;
    for (i = 0; i < N; i++)
    {
        
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Energy(M, dx, a2, a1, g, V, S);

        if ( i % 50 == 0 )
        {
            printf(" \n  %.4lf          ", i*dt);
            printf("%15.7E          ", creal(E));
            printf("%15.7E          ", Rsimps(M, abs2, dx));
        }



        // Apply exponential with nonlinear part
        rcarrExp(M, g * Idt / 2, abs2, exp_pot);
        carrMultiply(M, exp_pot, S, linpart);

        // Solve linear part
        CCSvec(m, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicLU(m, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        // Apply exponential with nonlinear part again
        carrAbs2(M, linpart, abs2);
        rcarrExp(M, g * Idt / 2, abs2, exp_pot);
        carrMultiply(M, exp_pot, linpart, S);


        // record data every n steps
        if (k == n) { carr_inline(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }

    }

    carrAbs2(M, S, abs2);
    E = Energy(M, dx, a2, a1, g, V, S);
    printf(" \n  %.4lf          ", N*dt);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));
    
    sepline();

    fclose(out_data);

    free(exp_pot);
    free(linpart);
    free(abs2);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    CCSFree(cnmat);
}





void SSCNSM(EqDataPkg EQ, int N, double dt, int cyclic, Carray S,
     char fname[], int n)
{

/** Evolve the wave-function given an initial condition in S
  * that is overwritten at each time-step.  The  results are
  * recorded in a file named 'fname' on every 'n' steps. Use
  * Crank-Nicolson discretization scheme to  compute  linear
  * part of potential and derivatives of the PDE.   'Cyclic'
  * is a boolean argument define  whether  the  boundary  is
  * zero or periodic on last position point.             **/



    unsigned int
        k,
        i,
        M,
        m,
        j;

    double
        a2,
        dx,
        g;

    double complex
        E,
        a1,
        Idt;

    Rarray
        V,
        abs2;

    Carray
        exp_pot,
        linpart,
        upper,
        lower,
        mid,
        rhs;

    CCSmat
        cnmat;

    FILE
        * out_data;



    M = EQ->Mpos; // grid size
    m = M - 1;    // grid size excluding the boundaries

    abs2 = rarrDef(M);    // squared modulus of wave function at grid points
    exp_pot = carrDef(M); // nonlinear potential part exponential in SS
    linpart = carrDef(M); // linear part solution by Finite Differences

    // tridiagonal system with additional  cyclic  terms depending on
    // the boundary conditions given.  In  positive case for periodic
    // boundary use the last element of 'upper' and 'lower' diagonals
    // to store the last column first line element in 'upper' and the
    // last line first column element in 'lower'.
    upper = carrDef(m);
    lower = carrDef(m);
    mid   = carrDef(m);

    rhs   = carrDef(m); // RHS of linear system to solve

    out_data = fopen(fname, "w");
    if (out_data == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    carr_inline(out_data, M, S);

    // unpack equation parameters from structure
    a2 = EQ->a2;
    a1 = EQ->a1;
    dx = EQ->dx;
    g = EQ->inter;
    V = EQ->V;
    Idt = 0.0 - dt * I;



            /*  ===========================================

                Setup RHS matrix from Crank-Nicolson Scheme

                ===========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(m, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(m, upper, dt / 2, V, mid);

    // fill upper diagonal
    carrFill(m, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[m-1] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[m-1] = 0;                                        }

    // fill lower diagonal
    carrFill(m, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[m-1] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[m-1] = 0;                                        }

    // Store in Compressed-Comlumn-Storage format
    cnmat = cyclic2CCS(m, upper, lower, mid);



          /*  ================================================

              Setup Cyclic tridiagonal matrix of linear system

              ================================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(m, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(m, upper, -dt / 2, V, mid);

    // fill upper diagonal
    carrFill(m, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[m-1] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[m-1] = 0;                                          }

    // fill lower diagonal
    carrFill(m, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[m-1] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[m-1] = 0;                                          }



    // Header of the screen output
    printf("\n\n\n");
    printf("     Time            Energy                   Norm");
    sepline();



    k = 1;
    for (i = 0; i < N; i++)
    {
        
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Energy(M, dx, a2, a1, g, V, S);

        if ( i % 50 == 0 )
        {
            printf(" \n  %.4lf          ", i*dt);
            printf("%15.7E          ", creal(E));
            printf("%15.7E          ", Rsimps(M, abs2, dx));
        }



        // Apply exponential with nonlinear part
        rcarrExp(M, g * Idt / 2, abs2, exp_pot);
        carrMultiply(M, exp_pot, S, linpart);

        // Solve linear part
        CCSvec(m, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicSM(m, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        // Apply exponential with nonlinear part again
        carrAbs2(M, linpart, abs2);
        rcarrExp(M, g * Idt / 2, abs2, exp_pot);
        carrMultiply(M, exp_pot, linpart, S);


        // record data every n steps
        if (k == n) { carr_inline(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }

    }

    carrAbs2(M, S, abs2);
    E = Energy(M, dx, a2, a1, g, V, S);
    printf(" \n  %.4lf          ", N*dt);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));
    
    sepline();

    fclose(out_data);

    free(exp_pot);
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

    Following the split step approach, the nonlinear part for the next
    methods is evolved in time using the 4th order Runge-Kutta  method
    since in a discretized grid it generates a system of ODEs. Our API
    of RK4 integrator for a system of complex ODE is rather general and
    can be find in a separate file, where we only need to provide here
    the derivative function required as paramenter there.

    =====================================================================  */





void NonLinearDDT(int M, double t, Carray Psi, Carray inter, Carray Dpsi)
{

/** In the output parameter 'Dpsi' compute derivative for equation given by
  * i df/dt = g|f(Xi,t)|^2 f(Xi,t)                                      **/

    int
        i;

    double
        mod_squared;

    for (i = 0; i < M; i++)
    {
        mod_squared  = creal(Psi[i]) * creal(Psi[i]);
        mod_squared += cimag(Psi[i]) * cimag(Psi[i]);
        Dpsi[i] = - I * inter[0] * mod_squared * Psi[i];
    }
}





void NonLinearVDDT(int M, double t, Carray Psi, Carray FullPot, Carray Dpsi)
{

/** In the output parameter 'Dpsi' compute derivative for equation given by
  * i df/dt = ( g|f(Xi,t)|^2 + V(Xi) ) f(Xi,t)                          **/

    int
        i;

    double
        mod_squared;

    for (i = 0; i < M; i++)
    {
        mod_squared  = creal(Psi[i]) * creal(Psi[i]);
        mod_squared += cimag(Psi[i]) * cimag(Psi[i]);
        Dpsi[i] = - I * (FullPot[0] * mod_squared + FullPot[i + 1] ) * Psi[i];
    }
}





void SSCNRK4(EqDataPkg EQ, int N, double dt, int cyclic, Carray S,
     char fname[], int n)
{

/** Similar to SSCN routine but use RK4 to evolve nonliear part **/



    int k,
        i,
        m,
        M,
        j;

    FILE
        * out_data;

    double
        a2,
        dx;

    double complex
        E,
        a1,
        g[1];

    Rarray
        V,
        abs2;

    Carray
        linpart,
        upper,
        lower,
        mid,
        rhs;

    CCSmat
        cnmat;



    M = EQ->Mpos;
    m = M - 1;

    a2 = EQ->a2;
    a1 = EQ->a1;
    dx = EQ->dx;
    g[0] = EQ->inter;
    V = EQ->V;

    abs2 = rarrDef(M);    // squared modulus of wave function at grid points
    linpart = carrDef(M); // linear part solution by Finite Differences

    // tridiagonal system with additional  cyclic  terms depending on
    // the boundary conditions given.  In  positive case for periodic
    // boundary use the last element of 'upper' and 'lower' diagonals
    // to store the last column first line element in 'upper' and the
    // last line first column element in 'lower'.
    upper = carrDef(m);
    lower = carrDef(m);
    mid   = carrDef(m);

    rhs   = carrDef(m); // RHS of linear system to solve

    out_data = fopen(fname, "w");
    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("\n\nERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    carr_inline(out_data, M, S);



    // Reader of screen printing
    printf("\n\n\n");
    printf("     time            Energy                   Norm");
    sepline();



            /*  ==========================================

                       Setup Right-Hand-Side matrix

                ==========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(m, - a2 * dt / dx / dx + I, upper);
    rcarrUpdate(m, upper, dt / 2, V, mid);

    // fill upper diagonal
    carrFill(m, a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, upper);
    if (cyclic) { upper[m-1] = a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { upper[m-1] = 0;                                        }

    // fill lower diagonal
    carrFill(m, a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, lower);
    if (cyclic) { lower[m-1] = a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { lower[m-1] = 0;                                        }

    // Store in CCS format
    cnmat = cyclic2CCS(m, upper, lower, mid);



            /*  ===========================================

                      Setup Cyclic tridiagonal matrix

                ===========================================  */



    // fill main diagonal (use upper as auxiliar pointer)
    carrFill(m, a2 * dt / dx /dx + I, upper);
    rcarrUpdate(m, upper, -dt / 2, V, mid);

    // fill upper diagonal
    carrFill(m, - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4, upper);
    if (cyclic) { upper[m-1] = - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4; }
    else        { upper[m-1] = 0;                                          }

    // fill lower diagonal
    carrFill(m, - a2 * dt / dx / dx / 2 + a1 * dt / dx / 4, lower);
    if (cyclic) { lower[m-1] = - a2 * dt / dx / dx / 2 - a1 * dt / dx / 4; }
    else        { lower[m-1] = 0;                                          }



    k = 1;
    for (i = 0; i < N; i++)
    {
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Energy(M, dx, a2, a1, g[0], V, S);
        if ( i % 50 == 0 )
        {
            printf(" \n  %.4lf          ", i*dt);
            printf("%15.7E          ", creal(E));
            printf("%15.7E          ", Rsimps(M, abs2, dx));
        }



        RK4step(M, dt/2, 0, S, g, linpart, NonLinearDDT);
        
        // Solve linear part (nabla ^ 2 part + onebody potential)
        CCSvec(m, cnmat->vec, cnmat->col, cnmat->m, linpart, rhs);
        triCyclicSM(m, upper, lower, mid, rhs, linpart);
        if (cyclic) { linpart[M-1] = linpart[0]; } // Cyclic system
        else        { linpart[M-1] = 0;          } // zero boundary

        RK4step(M, dt/2, 0, linpart, g, S, NonLinearDDT);



        // record data every n steps
        if (k == n) { carr_inline(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }
    }

    carrAbs2(M, S, abs2);
    E = Energy(M, dx, a2, a1, g[0], V, S);
    printf(" \n  %.4lf          ", N*dt);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));

    sepline();

    fclose(out_data);

    free(linpart);
    free(upper);
    free(lower);
    free(mid);
    free(rhs);
    free(abs2);
    CCSFree(cnmat);
}





void SSFFTRK4(EqDataPkg EQ, int N, double dt, Carray S, char fname[], int n)
{

/** Similar to SSFFT routine but uses RK4 to evolve nonlinear part **/



    int
        k,
        i,
        j,
        M,
        m;

    MKL_LONG
        s;

    double
        a2,
        dx,
        g,
        freq;

    double complex
        E,
        a1,
        Idt = 0 - dt * I;

    Rarray
        V,
        abs2;

    Carray
        argRK4,      // output of linear and nonlinear pot
        exp_der,     // exponential of derivative operator
        forward_fft, // go to frequency space
        back_fft,    // back to position space
        PotArg;

    FILE
        * out_data;

    DFTI_DESCRIPTOR_HANDLE
        desc;



    M = EQ->Mpos;
    m = EQ->Mpos - 1;



    out_data = fopen(fname, "w");
    if (out_data == NULL)
    {
        printf("\n\nERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    carr_inline(out_data, M, S);

    abs2 = rarrDef(M);

    argRK4 = carrDef(M);      // input in RK4
    exp_der = carrDef(m);     // exponential of derivative operator
    forward_fft = carrDef(m); // go to frequency space
    back_fft = carrDef(m);    // back to position space
    PotArg = carrDef(M + 1);  // Extra arguments for RK4



    a2 = EQ->a2;
    a1 = EQ->a1;
    dx = EQ->dx;
    g = EQ->inter;
    V = EQ->V;



    // Setup Extra arguments needed in RK4 routine
    PotArg[0] = g;
    for (i = 0; i < M; i++) PotArg[i + 1] = V[i];
    


    // setup descriptor (MKL implementation of FFT)
    s = DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, m);
    s = DftiSetValue(desc, DFTI_FORWARD_SCALE, 1.0 / sqrt(m));
    s = DftiSetValue(desc, DFTI_BACKWARD_SCALE, 1.0 / sqrt(m));
    s = DftiCommitDescriptor(desc);



    // setup Fourier Frequencies and the exponential of derivative operator
    for (i = 0; i < m; i++) {
        if (i <= (m - 1) / 2) { freq = (2 * PI * i) / (m * dx);       }
        else                  { freq = (2 * PI * (i - m)) / (m * dx); }
        // exponential of derivative operators
        exp_der[i] = cexp(Idt * a1 * freq * I - Idt * a2 * freq * freq);
    }



    // Header of screen output
    printf("\n\n\n");
    printf("     time            Energy                   Norm");
    sepline();



    k = 1;
    for (i = 0; i < N; i++)
    {

        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Energy(M, dx, a2, a1, g, V, S);
        if ( i % 50 == 0 )
        {
            printf(" \n  %.4lf          ", i*dt);
            printf("%15.7E          ", creal(E));
            printf("%15.7E          ", Rsimps(M, abs2, dx));
        }



        RK4step(M, dt/2, 0, S, PotArg, argRK4, NonLinearVDDT);
        carrCopy(m, argRK4, forward_fft);

        // go to momentum space
        s = DftiComputeForward(desc, forward_fft);
        // apply exponential of derivatives
        carrMultiply(m, exp_der, forward_fft, back_fft);
        // go back to position space
        s = DftiComputeBackward(desc, back_fft);
        carrCopy(m, back_fft, argRK4);
        argRK4[m] = argRK4[0]; // cyclic condition

        RK4step(M, dt/2, 0, argRK4, PotArg, S, NonLinearVDDT);

        // record data every n steps
        if (k == n) { carr_inline(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }

    }

    carrAbs2(M, S, abs2);
    E = Energy(M, dx, a2, a1, g, V, S);
    printf(" \n  %.4lf          ", N*dt);
    printf("%15.7E          ", creal(E));
    printf("%15.7E          ", Rsimps(M, abs2, dx));

    sepline();

    fclose(out_data);

    s = DftiFreeDescriptor(&desc);

    free(exp_der);
    free(forward_fft);
    free(back_fft);
    free(argRK4);
    free(abs2);
    free(PotArg);
}










void CFDS(EqDataPkg EQ, int N, double dt, int cyclic, Carray S,
     char fname [], int n)
{



    unsigned int
        k,
        i,
        j,
        M,
        iter,
        condition;

    M = EQ->Mpos;


    // File to write every n step the time-step solution
    FILE * out_data = fopen(fname, "w");



    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("\n\n\tERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    carr_inline(out_data, M, S);



    // Reader of screen printing
    printf("\n\n\n");
    printf("     step            Energy                   Norm");
    sepline();



    double
        a2,
        dx,
        inter,
        * V,
        tol;



    // Factor to multiply in nonlinear exponential part after split-step
    double complex
        a1,
        aux,
        Idt = 0.0 - dt * I;



    Rarray
        abs2 = rarrDef(M); // abs square of wave function



    Carray
        // Fixed point iteratively solve
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



    a2 = EQ->a2;
    a1 = EQ->a1;
    dx = EQ->dx;
    inter = EQ->inter;
    V = EQ->V;



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
    cnmat = cyclic2CCS(M - 1, upper, lower, mid);



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
        aux = Energy(M, dx, a2, a1, inter, V, S);
        if ( i % 50 == 0 )
        {
            printf(" \n  %7d          ", i);
            printf("%15.7E          ", creal(aux));
            printf("%15.7E          ", Rsimps(M, abs2, dx));
        }



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
                tol = 1E-13 * cabs(linpart[j]) + 1E-13;
                if (cabs(Sstep[j] - linpart[j]) > tol) condition = 1;
            }

            iter = iter + 1;

        }

        // correct the main diagonal to do SSCN
        for (j = 0; j < M - 1; j ++)
        {
            mid[j] = I + a2 * dt / dx / dx - dt * V[j] / 2;
        }

        carrCopy(M, Sstep, S);
        
        // record data every n steps
        if (k == n) { carr_inline(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }

    }

    carrAbs2(M, S, abs2);
    aux = Energy(M, dx, a2, a1, inter, V, S);
    printf(" \n  %7d          ", N);
    printf("%15.7E          ", creal(aux));
    printf("%15.7E          ", Rsimps(M, abs2, dx));

    sepline();

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





void sinedvrDDT(int M, double t, Carray a, Carray args, Carray Da)
{

/** In the output parameter 'Dpsi' compute derivative for equation given by
  * i df/dt = ( g|f(Xi,t)|^2 + V(Xi) ) f(Xi,t)                          **/

    int
        i,
        j;

    double
        w,
        g;

    double complex
        sum;

    w = creal(args[M * M]);     // weights of quadrature
    g = creal(args[M * M + 1]); // interaction strength

    for (i = 0; i < M; i++)
    {
        sum = 0;
        for (j = 0; j < M; j++)
        {
            sum = sum + args[i*M + j] * a[j];
        }

        sum = sum + g * conj(a[i]) * a[i] * a[i] / w;

        Da[i] = - I * sum;
    }
}





void sineDVR(EqDataPkg EQ, int N, double dt, Carray S, char fname[], int n)
{



    int
        l,
        k,
        i,
        M,
        j,
        i1,
        j1;

    FILE
        * out_data;

    double
        L,
        g,
        a2,
        dx,
        sum;

    double complex
        E,
        a1;

    Rarray
        abs2,
        D1mat,
        D2mat,
        uDVR;

    Carray
        aDVR,
        D2DVR,
        D1DVR,
        RK4arg;



    M = EQ->Mpos;        // number of functions in the basis and consequently
                         // the number of points to represent the solution

    // unpack equation coefficients and array with potential at grid points
    // In sine DVR the grid step is the weights of gaussian quadrature
    a2 = EQ->a2;
    a1 = EQ->a1;
    dx = EQ->dx;
    g = EQ->inter;

    L = (M + 1) * dx; // length of the sine DVR box

    abs2 = rarrDef(M);

    D1mat = rarrDef(M * M); // derivative matrix
    D2mat = rarrDef(M);     // second derivative matrix(diagonal)

    D1DVR = carrDef(M * M); // derivative matrix in DVR basis
    D2DVR = carrDef(M * M); // second derivative matrix in DVR basis

    uDVR = rarrDef(M * M);  // unitary transformation to DVR basis

    aDVR = carrDef(M);  // wave function coefficients in DVR basis

    RK4arg = carrDef(M * M + 2);



    out_data = fopen(fname, "w");
    if (out_data == NULL)
    {   // impossible to open file with the given name
        printf("\n\nERROR: impossible to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }

    // Record initial data as first line
    carr_inline(out_data, M, S);



    // Setup DVR transformation matrix with  eigenvector  organized
    // by columns in row major format. It also setup the derivarive
    // matrices in the sine basis that are known analytically
    for (i = 0; i < M; i++)
    {
        i1 = i + 1;
        for (j = 0; j < M; j++)
        {
            j1 = j + 1;
            uDVR[i*M + j] = sqrt(2.0 / (M + 1)) * sin(i1*j1*PI/(M+1));
            if ( (i1 - j1) % 2 != 0)
            {
                D1mat[i*M + j] = 2.0 * j1 / L / (i1 - j1);
            }
            else
            {
                D1mat[i*M + j] = 0;
            }
        }
        D2mat[i] = - (i1 * PI / L) * (i1 * PI / L);
    }



    printf("\n\nTransform Matrix to DVR basis\n");

    if (cabs(a1) != 0)
    {
        // Transform first order derivative matrix to DVR basis
        // multiplying by the equation coefficients
        for (i = 0; i < M; i++)
        {
            printf("\nDoing line %d",i);
            for (j = 0; j < M; j++)
            {
                sum = 0;
                for (k = 0; k < M; k++)
                {
                    for (l = 0; l < M; l++)
                    {
                        sum += uDVR[i*M + k] * D1mat[k*M + l] * uDVR[l*M + j];
                    }
                }
                D1DVR[i*M + j] = a1 * sum;
            }
        }
    }
    else
    {
        printf("\nNo need of first order derivative\n");
        for (i = 0; i < M; i++)
        {
            for (j = 0; j < M; j++) D1DVR[i*M + j] = 0;
        }
    }



    // Transform second order derivative matrix to DVR basis
    // multiplying by the equation coefficient
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < M; j++)
        {
            sum = 0;
            for (k = 0; k < M; k++)
            {
                sum += uDVR[i*M + k] * D2mat[k] * uDVR[k*M + j];
            }
            D2DVR[i*M + j] = a2 * sum;
        }
    }

    printf("\n\nDVR successfully assembled\n\n");



    // configure extra parameters for RK4 routine, those that are
    // constant in time evolution, all terms but nonlinear part
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < M; j++)
        {
            RK4arg[i*M + j] = D1DVR[i*M + j] + D2DVR[i*M + j];
        }
        RK4arg[i*M + i] = RK4arg[i*M + i] + EQ->V[i];
    }

    RK4arg[M * M] = dx;    // the weight is required for nonlinear part
    RK4arg[M * M + 1] = g; // contact interaction strength



    // setup the DVR coefficients of the wave-function from basis expansion
    // In this case the weights are constants 'dx'
    for (i = 0; i < M; i++)
    {
        aDVR[i] = sqrt(dx) * S[i];
    }



    // For the following operations, the matrices of derivatives in both
    // basis as well as the basis transformation are no longer required
    free(D2mat);
    free(D1mat);
    free(D1DVR);
    free(D2DVR);
    free(uDVR);



    // Header of screen printing
    printf("\n\n\n");
    printf("     time            Energy                   Norm");
    sepline();



    k = 1;
    for (i = 0; i < N; i++)
    {
        carrAbs2(M, S, abs2);

        // Print in screen to quality and progress control
        E = Energy(M, dx, a2, a1, g, EQ->V, S);
        if ( i % 50 == 0 )
        {
            printf(" \n  %.4lf          ", i*dt);
            printf("%15.7E          ", creal(E));
            printf("%15.7E          ", Rsimps(M,abs2,dx));
        }



        RK4step(M,dt,i*dt,aDVR,RK4arg,S,sinedvrDDT);
        // Transform back the solution to grid points
        // and update initial condition to next step
        carrCopy(M,S,aDVR);
        for( j = 0; j < M; j++)
        {
            S[j] = S[j] / sqrt(dx);
        }



        // record data every n steps
        if (k == n) { carr_inline(out_data, M, S); k = 1; }
        else        { k = k + 1;                          }
    }

    carrAbs2(M, S, abs2);
    E = Energy(M, dx, a2, a1, g, EQ->V, S);
    printf(" \n  %.4lf          ",N*dt);
    printf("%15.7E          ",creal(E));
    printf("%15.7E          ",Rsimps(M,abs2,dx));

    sepline();

    fclose(out_data);

    free(abs2);
    free(aDVR);
    free(RK4arg);
}
