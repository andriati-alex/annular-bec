#include "../include/NewtonCG.h"

void iteration_info(struct IterNCG It)
{
    printf("\n\n\t\t*** Analysis of the number of Iterations ***\n");

    printf("\n\t%d Newton iterations", It.newton);
    printf("\n\t%d Max CG iterations in a Newton loop", It.cg_max);
    printf("\n\t%d Min CG iterations in a Newton loop", It.cg_min);
    printf("\n\t%d CG iterations altogether", It.cg_total);
    printf("\n\tMean of %.1f CG iterations per Newton loop", It.cg_mean);
    printf("\n\tStd  +- %.1f in CG iterations\n", It.cg_std);
}

struct IterNCG ncg(int M, double tol, int maxiter, double dx, double a2,
                   double complex a1, double inter, double mu, Rarray V,
                   Carray f0)
{
    /*** Structure to analyse convergence ***/

    struct IterNCG iterations;
    iterations.cg_total = 0;
    iterations.cg_mean = 0;
    iterations.cg_max = 0;
    iterations.cg_min = INT_MAX;
    iterations.cg_std = 0;

    int * cgiterations = (int * ) malloc(maxiter * sizeof(int));

    int j,
        ni,         // Number of Newton iterations
        cgi,        // Count number of CG iterations needed in a Newton loop
        N = 2 * M;  // N the size of the concatenated system (real + imag)

    double EPS, eps;  // Accepted Error for Newton and CG methods respectively

    /***     First order derivative with purely imaginary coeficient     ***/

    double ddx = cimag(a1) / (2 * dx);

    /***           Terms that come from second derivative part           ***/

    double diag = - 2 * a2 / (dx * dx), // main diagonal values
           offd = a2 / (dx * dx);       // off  diagonal values

    /***              Terms that come from interaction part              ***/

    double int_rr, // variation of real part of real equation
           int_ri, //     ''    of imag part of real equation
           int_ii; //     ''    of imag part of imag equation

    RCCSmat A = emptyRCCS(N, 6); // Six non-zero elements per line

    Carray L0f = carrDef(M);  // Complex Right hand side
    Rarray rhs = rarrDef(N);  // Right hand side of linear system
    Rarray abs2 = rarrDef(M);

    Rarray xcg = rarrDef(N); // Solution of linear system by CG method

    /*** Tri-diagonal system to use as preconditioning ***/

    Rarray upper = rarrDef(N - 1);
    Rarray lower = rarrDef(N - 1);
    Rarray mid = rarrDef(N);

    rarrFill(N - 1, offd, upper);
    rarrFill(N - 1, offd, lower);
    rarrFill(N, diag - mu, mid);
    upper[M - 1] = - ddx;
    lower[M - 1] = - ddx;

    /*      ********************************************************      */
    /*       CONSTANT PART OF THE SYSTEM THAT ARISES FROM REAL PART       */
    /*      ********************************************************      */

    /*** setup first line ***/

    setValueRCCS(N, 0, 1, 1, offd, A);
    setValueRCCS(N, 0, 2, M - 1, offd, A);
    // First order derivative couples real to imaginary system
    setValueRCCS(N, 0, 4, M + 1, - ddx, A);
    setValueRCCS(N, 0, 5, N - 1, ddx, A);

    /*** setup last line ***/

    setValueRCCS(N, M - 1, 0, 0, offd, A);
    setValueRCCS(N, M - 1, 1, M - 2, offd, A);
    // First order derivative couples real to imaginary
    setValueRCCS(N, M - 1, 3, M, - ddx, A);
    setValueRCCS(N, M - 1, 4, N - 2, ddx, A);

    /** Kernel (no cyclic terms) **/
    for (j = 1; j < M - 1; j++) {
        setValueRCCS(N, j, 0, j - 1, offd, A);
        setValueRCCS(N, j, 2, j + 1, offd, A);
        // First Order derivative
        setValueRCCS(N, j, 3, j - 1 + M, ddx, A);
        setValueRCCS(N, j, 5, j + 1 + M, - ddx, A);
    }

    /*      *********************************************************      */
    /*     CONSTANT PART OF THE SYSTEM THAT ARISES FROM IMAGINARY PART     */
    /*      *********************************************************      */

    /*** setup first line ***/

    setValueRCCS(N, M, 4, M + 1, offd, A);
    setValueRCCS(N, M, 5, N - 1, offd, A);
    // First order derivative couples real to imaginary
    setValueRCCS(N, M, 1, 1, ddx, A);
    setValueRCCS(N, M, 2, M - 1, - ddx, A);

    /*** setup last line ***/

    setValueRCCS(N, N - 1, 3, M, offd, A);
    setValueRCCS(N, N - 1, 4, N - 2, offd, A);
    // First order derivative couples real to imaginary
    setValueRCCS(N, N - 1, 0, 0, ddx, A);
    setValueRCCS(N, N - 1, 1, M - 2, - ddx, A);

    /** Kernel (no cyclic terms) **/
    for (j = M + 1; j < N - 1; j++) {
        setValueRCCS(N, j, 3, j - 1, offd, A);
        setValueRCCS(N, j, 5, j + 1, offd, A);
        // First Order derivative
        setValueRCCS(N, j, 0, j - 1 - M, - ddx, A);
        setValueRCCS(N, j, 2, j + 1 - M, ddx, A);
    }

    /*     *********************************************************     */
    /*        ITERATION DEPEDENT PART ARISES FROM NON-LINEARITIES        */
    /*     *********************************************************     */

    applyL0(M, f0, dx, a2, a1, V, inter, mu, L0f);
    EPS = carrMod(M, L0f);
    
    printf("\nModulus of L0f = Sqrt(sum of residues squared):\n");

    printf("\n\t%.2E in Initial Guess\n", EPS);

    ni = 1; // start Newton iteration counter
    while (EPS > tol)
    {
        if (ni > maxiter) {
            printf("\n\n\tEXCEEDED MAXIMUM NUMBER OF ITERATIONS.\n\n");
            exit(-1);
        }

        carrAbs2(M, f0, abs2);
        rarrFill(N, 0, xcg); // Aways starts with zero as a guess for CG

        /* Add-up diagonal terms to contribution of second order derivative */

        // First Line
        int_rr = inter * ( abs2[0] + 2 * creal(f0[0]) * creal(f0[0]) );
        int_ii = inter * ( abs2[0] + 2 * cimag(f0[0]) * cimag(f0[0]) );
        setValueRCCS(N, 0, 0, 0, diag + V[0] - mu + int_rr, A);
        setValueRCCS(N, M, 3, M, diag + V[0] - mu + int_ii, A);

        // Last Line
        int_rr = inter * ( abs2[M-1] + 2 * creal(f0[M-1]) * creal(f0[M-1]) );
        int_ii = inter * ( abs2[M-1] + 2 * cimag(f0[M-1]) * cimag(f0[M-1]) );
        setValueRCCS(N, M - 1, 2, M - 1, diag + V[M-1] - mu + int_rr, A);
        setValueRCCS(N, N - 1, 5, N - 1, diag + V[M-1] - mu + int_ii, A);

        // Kernel (no cyclic terms)
        #pragma omp parallel for private(j, int_rr, int_ii, int_ri)
        for (j = 1; j < M - 1; j++) {
            int_rr = inter * ( abs2[j] + 2 * creal(f0[j]) * creal(f0[j]) );
            setValueRCCS(N, j, 1, j, diag + V[j] - mu + int_rr, A);
            int_ii = inter * ( abs2[j] + 2 * cimag(f0[j]) * cimag(f0[j]) );
            setValueRCCS(N, j + M, 4, j + M, diag + V[j] - mu + int_ii, A);
            int_ri = 2 * inter * cimag(f0[j]) * creal(f0[j]);
            setValueRCCS(N, j, 4, j + M, int_ri, A); // Coupling real to imag
            setValueRCCS(N, j + M, 1, j, int_ri, A); // Coupling imag to real
        }

        int_ri = 2 * inter * cimag(f0[0]) * creal(f0[0]);
        setValueRCCS(N, 0, 3, M, int_ri, A); // Coupling real to imag
        setValueRCCS(N, M, 0, 0, int_ri, A); // Coupling imag to real

        int_ri = 2 * inter * cimag(f0[M - 1]) * creal(f0[M - 1]);
        setValueRCCS(N, M - 1, 5, N - 1, int_ri, A); // Coupling real to imag
        setValueRCCS(N, N - 1, 2, M - 1, int_ri, A); // Coupling imag to real

        /*** End of matrix entries ***/

        for (j = 0; j < M; j++) {
            rhs[j] = - creal(L0f[j]);
            rhs[j + M] = - cimag(L0f[j]);
        }

        eps = 0.005 * EPS;
        cgi = RCG(N, A, rhs, xcg, eps, upper, lower, mid);
        cgiterations[ni - 1] = cgi;

        /*** update Solution ***/

        for (j = 0; j < M; j++) f0[j] += xcg[j] + I * xcg[j + M];

        applyL0(M, f0, dx, a2, a1, V, inter, mu, L0f);
        EPS = carrMod(M, L0f);

        printf("\n\t%.2E in Newton iteration #%d", EPS, ni);
        printf("\n\t%d CG iterations\n", cgi);

        if (cgi > iterations.cg_max) iterations.cg_max = cgi;
        if (cgi < iterations.cg_min) iterations.cg_min = cgi;
        iterations.cg_mean  += cgi;
        iterations.cg_total += cgi;

        ni += 1; // Update Newton iteration counter
    }

    iterations.cg_mean /= (ni - 1);

    for (j = 0; j < ni - 1; j++)
        iterations.cg_std += pow((iterations.cg_mean - cgiterations[j]), 2);
    iterations.cg_std = sqrt(iterations.cg_std / (ni - 1));

    iterations.newton = ni;

    /*** Release used memory ***/

    free(abs2); RCCSFree(A); free(upper); free(lower); free(mid); free(xcg);
    free(rhs); free(L0f); free(cgiterations);

    return iterations;
}
