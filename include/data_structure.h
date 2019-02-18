#ifndef _data_structure_h
#define _data_structure_h

#include <string.h>
#include "array_memory.h"
#include "linear_potential.h"



/* ======================================================================== */
/*                                                                          */
/*         DATA-TYPE DEFINITION - All information Needed for MCTDHB         */
/*                                                                          */
/* ======================================================================== */



struct _EquationDataPkg
{

    char
        Vname[80]; // One-Body potential

    int 
        Mpos;     // # of discretized positions (# divisions + 1)

    double
        dx,       // space step
        xi,       // initial position discretized value
        xf,       // final position discretized value
        a2,       // factor multiplying d2 / dx2
        inter,    // know as g, contact interaction strength
        * V;      // Array with the values of one-particle potential

    double
        p[3];     // Parameters to generate one-particle potential values

    double complex
        a1;       // factor multiplying d / dx (pure imaginary)

};

typedef struct _EquationDataPkg * EqDataPkg;



/* ======================================================================== *
 *                                                                          *
 *                           FUNCTION PROTOTYPES                            *
 *                                                                          *
 * ======================================================================== */





EqDataPkg PackEqData(int,double,double,double,double,doublec,char [],double []);

void ReleaseEqDataPkg (EqDataPkg);



#endif
