#ifndef _MCTDHB_datatype_h
#define _MCTDHB_datatype_h

#include <complex.h>
#include "MCTDHB_configurations.h"



/* ======================================================================== */
/*                                                                          */
/*         DATA-TYPE DEFINITION - All information Needed for MCTDHB         */
/*                                                                          */
/* ======================================================================== */



struct _MCTDHBsetup
{
    int 
        Mpos,   // # of discretized positions (# divisions + 1)
        Morb,   // # of orbitals
        Npar,   // # of particles
        ** IF;  // IF[i] point to the occupation number vetor of C[i]
    int
        nc,         // Total # of configurations of Fock states
        ** NCmat;   // NCmat[n][m] # with n particles / m orbitals
    double
        dx,     // space step
        xi,     // initial position discretized value
        xf,     // final position discretized value
        a2,     // factor multiplying d2 / dx2
        inter,  // know as g, contact interaction strength
        * V;    // Array with the values of one-particle potential
    double complex
        a1;     // factor multiplying d / dx (pure imaginary)
};

typedef struct _MCTDHBsetup * MCTDHBsetup;



/* ======================================================================== *
 *                                                                          *
 *                           FUNCTION PROTOTYPES                            *
 *                                                                          *
 * ======================================================================== */



MCTDHBsetup AllocMCTDHBdata
(
    int Npar,
    int Morb,
    int Mpos,
    double xi,
    double xf,
    double a2,
    double inter,
    double * V,
    double complex a1
);

void EraseMCTDHBdata(MCTDHBsetup MC);

#endif
