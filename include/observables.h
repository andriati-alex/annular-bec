#ifndef _observables_h
#define _observables_h

#include "calculus.h"





doublec KinectE(int M, double a2, doublec a1, double dx, Carray psi);



doublec TrapE(int M, Rarray V, double dx, Carray psi);



double InterE(int M, double g, double dx, Carray psi);



doublec Energy(int M, double dx, double a2, doublec a1, double inter,
        Rarray V, Carray f);



doublec Chem(int M, double dx, double a2, doublec a1, double inter,
        Rarray V, Carray f);



doublec Virial(int M, double a2, doublec a1, double g, Rarray V,
        double dx, Carray psi);





/** For distributions in symmetric domain **/
double MeanQuadraticR(int n, Carray f, double dx);





#endif
