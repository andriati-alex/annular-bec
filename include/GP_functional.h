#ifndef _GP_functional_h
#define _GP_functional_h

#include "calculus.h"





double complex Functional(int M, double dx, double a2, double complex a1,
               double inter, Rarray V, Carray f);



void applyL0(int n, Carray f, double dx, double a2, double complex a1, 
             Rarray V, double inter, double mu, Carray L0f);



double complex GPkinect(int M, double a2, double complex a1, double dx,
               Carray psi);



double complex GPtrap(int M, Rarray V, double dx, Carray psi);



double GPinter(int M, double g, double dx, Carray psi);



double complex GPvirial(int M, double a2, double complex a1, double g,
               Rarray V, double dx, Carray psi);





#endif
