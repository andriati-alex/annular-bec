#ifndef _observables_h
#define _observables_h

#include "calculus.h"





doublec Functional(int M, double dx, double a2, doublec a1, double inter,
        Rarray V, Carray f);



void applyL0(int n, Carray f, double dx, double a2, double complex a1, 
             Rarray V, double inter, double mu, Carray L0f);



doublec kinectE(int M, double a2, doublec a1, double dx, Carray psi);



doublec TrapE(int M, Rarray V, double dx, Carray psi);



double InterE(int M, double g, double dx, Carray psi);



doublec Virial(int M, double a2, doublec a1, double g, Rarray V,
        double dx, Carray psi);





#endif
