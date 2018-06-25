#include "../include/NewtonCG.h"

void ncg(int M, double dx, double, a2, double complex a1, double inter,
         Rarray V, Carray f0, Carray S, int * ni, int * cgi)
{

    int j, N = 2 * M; // N the size of the concatenated system

    /* Terms that come from derivative part */
    double diag   = - 2 * a2 / (dx * dx), // main diagonal values
           offd = a2 / (dx * dx);         // off diagonal values

    /* Terms that come from interaction part */
    double int_rr, // variation of real part of real equation
           int_ri, //     ''    of imag part of real equation
           int_ii; //     ''    of imag part of imag equation

    RCCSmat A = emptyCCS(N, 4);

    Rarray abs2 = rarrDef(M);

    /*     *********************************************************     */
    /*        FIRST PART OF THE SYSTEM THAT ARISES FROM REAL PART        */
    /*     *********************************************************     */

    /** setup first line **/
    A->vec[0] = diag;
    A->col[0] = 0;
    A->vec[0 + 1 * N] = offd;
    A->col[0 + 1 * N] = 1;
    A->vec[0 + 2 * N] = offd;
    A->col[0 + 2 * N] = M - 1;
    /** setup last line **/
    A->vec[M - 1] = offg;
    A->col[M - 1] = 0;
    A->vec[M - 1 + 1 * N] = offd;
    A->col[M - 1 + 1 * N] = M - 2;
    A->vec[M - 1 + 2 * N] = diag;
    A->col[M - 1 + 2 * N] = M - 1;
    /** Kernel (no cyclic terms) **/
    for (i = 1; i < M - 1; i++) {
        setValueRCCS(N, j, 1, j, diag, A);
        setValueRCCS(N, j, 0, j - 1, offd, A);
        setValueRCCS(N, j, 2, j + 1, offd, A);
    }

    /*     *********************************************************     */
    /*     SECOND PART OF THE SYSTEM THAT ARISES FROM IMAGINARY PART     */
    /*     *********************************************************     */

    /** setup first line **/
    A->vec[M + 1 * N] = diag;
    A->col[M + 1 * N] = M;
    A->vec[M + 2 * N] = offd;
    A->col[M + 2 * N] = M + 1;
    A->vec[M + 3 * N] = offd;
    A->col[M + 3 * N] = N;
    /** setup last line **/
    A->vec[N - 1 + 1 * N] = offg;
    A->col[N - 1 + 1 * N] = M;
    A->vec[N - 1 + 2 * N] = offd;
    A->col[N - 1 + 2 * N] = N - 2;
    A->vec[N - 1 + 3 * N] = diag;
    A->col[N - 1 + 3 * N] = N - 1;
    /** Kernel (no cyclic terms) **/
    for (j = M + 1; j < N - 1; j++) {
        setValueRCCS(N, j, 2, j, diag, A);
        setValueRCCS(N, j, 1, j - 1, offd, A);
        setValueRCCS(N, j, 3, j + 1, offd, A);
    }
    
    /*     *********************************************************     */
    /*        ITERATION DEPEDENT PART ARISES FROM NON-LINEARITIES        */
    /*     *********************************************************     */

    carrAbs2(M, f0, abs2);

    // Add-up terms to main diagonal
    int_rr = inter * ( abs2[0] + 2 * creal(f0[0]) * creal(f0[0]) );
    int_ii = inter * ( abs2[0] + 2 * cimag(f0[0]) * cimag(f0[0]) );
    addValueRCCS(N, 0, 0, int_rr, A);
    addValueRCCS(N, M, 1, int_ii, A);

    int_rr = inter * ( abs2[M-1] + 2 * creal(f0[M-1]) * creal(f0[M-1]) );
    int_ii = inter * ( abs2[M-1] + 2 * cimag(f0[M-1]) * cimag(f0[M-1]) );
    addValueRCCS(N, M - 1, 2, int_rr, A);
    addValueRCCS(N, N - 1, 3, int_ii, A);

    #pragma omp parallel for private(i)
    for (j = 1; j < M - 1; j++) {
        int_rr = inter * ( abs2[j] + 2 * creal(f0[j]) * creal(f0[j]) );
        addValueRCCS(N, j, 1, int_rr, A);
        int_ii = inter * ( abs2[j] + 2 * cimag(f0[j]) * cimag(f0[j]) );
        addValeuRCCS(N, j + M, 2, int_ii, A);
        int_ri = 2 * inter * cimag(f0[j]) * creal(f0[j]);
        setValueRCCS(N, j, 3, j + M, int_ri, A); // Coupling real to imag
        setValueRCCS(N, j + M, 0, j, int_ri, A); // Coupling imag to real
    }
    
    int_ri = 2 * inter * cimag(f0[0]) * creal(f0[0]);
    setValueRCCS(N, 0, 3, M, int_ri, A); // Coupling real to imag
    setValueRCCS(N, M, 0, 0, int_ri, A); // Coupling imag to real

    int_ri = 2 * inter * cimag(f0[M - 1]) * creal(f0[M - 1]);
    setValueRCCS(N, M - 1, 3, N, int_ri, A); // Coupling real to imag
    setValueRCCS(N, N, 0, M - 1, int_ri, A); // Coupling imag to real

    free(abs2); RCCSFree(A);
}
