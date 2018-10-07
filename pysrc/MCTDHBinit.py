
import sys;
import numpy as np;
import MCTDHBmodule as mc;

from scipy.integrate import simps;
from math import pi;
from math import sqrt;





"""
=============================================================================


    SCRIPT TO GENERATE INITIAL CONDITION DATA FOR MCTDHB
    ----------------------------------------------------

    Generate 3 files in setup/ folder:



    1) 'MC_fileId_config.dat' contains the numbers:

        Npar Morb Mdiv xi xf



    2) 'MC_fileId_orb.dat' contains:

        Matrix whose each column is an orbital and the values of
        the orbitals in discretized positions is given over  the
        rows. For example if Morb = 3 then:

                Orbital #1     Orbital #2     Orbital #3

           x1 |   f1(x1)         f2(x1)         f3(x1)
           x2 |   f1(x2)         f2(x2)         f3(x2)
              |     .              .              .
              |     .              .              .
              |     .              .              .
           xm |   f1(xm)         f2(xm)         f3(xm)



    3) 'MC_fileId_coef.dat' contains a column vector of size:

         (Npar + Morb - 1)!   | with values of coefficietns of
        --------------------  | all possible configurations on
        (Npar)!  (Morb - 1)!  | occupations in a Fock state


=============================================================================
"""





lf = np.float128;
lc = np.complex256;





def renormalize(f, dx): return f / sqrt( simps(abs(f)**2, dx = dx) );





def BrightSoliton(x, a, c):
    """ Normalized bright soliton with sech shape  """

    numerator = a * np.exp(0.5j * c * x * sqrt(2), dtype=lc);
    denominator = np.cosh(a * x / sqrt(2), dtype=lf);
    return numerator / denominator / np.sqrt(2 * sqrt(2) * a);





def AngularMom(Morb, x, S):
    """
    Setup matrix S with eigenstates of  angular  momentum  in ascending
    if the # of orbitals(Morb) is odd. Odd orbitals have minus signs of
    angular momentum while even orbitals have plus sign.
    -------------------------------------------------------------------

    x    = array of discretized position domain
    Morb = # of orbitals
    """

    S[0,:] = 1.0 / sqrt(2 * pi);
    for i in range(2, Morb, 2):
        S[i - 1,:] = np.exp(- 1.0j * (i / 2) * x, dtype=lc) / sqrt(2 * pi);
        S[i, :]    = np.exp(  1.0j * (i / 2) * x, dtype=lc) / sqrt(2 * pi);





def NoiseAngular(Morb, x, S, k):
    """
    Setup matrix 'S' with superposition of 'k' pairs of degenerated angular
    momentum eigenstates, i.e, the plus with the minus sign,  using  random
    number for the coefficient to do the combination.
    -----------------------------------------------------------------------

    x    = array of discretized position domain
    Morb = # of orbitals
    """

    phase = 0 + 0j;
    # n is used to index the random numbers
    n = 0;
    # Random phase and amplitude
    c = (np.random.random(2 * k * Morb) - 0.5) / 0.1;
    for i in range(Morb):
        for l in range(k):
            # add k pairs of angular momentum
            n = 2 * l + 2 * i * k;
            # positive angular quantum number
            phase = np.exp(- 1.0j * ( (i*k + l) * x + c[n]), dtype=lc);
            S[i,:] += c[n] * phase;
            # negative angular quantum number
            phase = np.exp(+ 1.0j * ( (i*k + l) * x + c[n + 1]), dtype=lc);
            S[i,:] += c[n + 1] * phase;
        S[i,:] = renormalize(S[i,:], x[1] - x[0]);





def ThermalCoef(Npar, Morb, beta, C):
    """
    Consider the energy proportional to the square  of  orbital  number
    and them consider a thermal-like distribution  with  the fock-state
    coefficient decreasing exponentially according  to  the  number  of
    occupation in each orbital(j) times the square of orbital number(j)
    (representing the energy)
    -------------------------------------------------------------------

    C = array of coefficient of size of # of all possible fock config.
    Npar = # of particles
    Morb = # of orbitals
    beta = decay rate (something like ~ 1 / temperature)
    """
    nc = mc.NC(Npar, Morb);
    v = np.empty(Morb, dtype=np.int32); # Fock vector for each coefficient
    phase = np.exp(2 * pi * 1.0j * (np.random.random(nc) - 0.5));

    for l in range(nc):
        mc.IndexToFock(l, Npar, Morb, v);
        prod = 1.0;
        for j in range(Morb):
            decay = - beta * float(v[j] * j * j) / Npar;
            prod  = prod * np.exp( decay , dtype=lf );
        # put the random phase with 'thermal' exponential decay
        C[l] = phase[l] * prod;

    # renormalize coefficients
    Norm = 0;
    for l in range(nc): Norm = Norm + abs(C[l]) * abs(C[l]);
    for l in range(nc): C[l] = C[l] / sqrt(Norm);





def BECCoef(Npar, Morb, C):
    """ All particles occupying the orbital 0 """

    nc = mc.NC(Npar, Morb);
    C[0] = 1.0;
    for i in range(1, nc): C[i] = 0;










###################      Read command line options      ####################

Npar = int(sys.argv[1]);   # Number of Particles
Morb = int(sys.argv[2]);   # Number of orbitals
xi   = float(sys.argv[3]); # initial position
xf   = float(sys.argv[4]); # final position
Mdiv = int(sys.argv[5]);   # Number of divisions slices in domain
Id   = int(sys.argv[6]);   # function to generate Id

params = []; # extra parameters to pass to the function
for i in range(7, len(sys.argv)): params.append( lf(sys.argv[i]) );
params = tuple(params);





############### Call subroutines to setup initial condition ################

x  = np.linspace(xi, xf, Mdiv + 1, dtype=lf);
dx = (xf - xi) / Mdiv;

Orb = np.zeros([Morb, x.size], dtype=lc);  # orbitals
C = np.zeros(mc.NC(Npar, Morb), dtype=lc); # coeficients

if (Id == 1) :
    Id_name = 'AngularMom-' + str(Npar) + '-' + str(Morb);
    AngularMom(Morb, x, Orb);
    ThermalCoef(Npar, Morb, params[0], C);
else :
    Id_name = 'NoiseAngular-' + str(Npar) + '-' + str(Morb);
    NoiseAngular(Morb, x, Orb, int(params[1]));
    ThermalCoef(Npar, Morb, params[0], C);





#################              Record Data              ###################

np.savetxt('../setup/MC_' + Id_name + '_orb.dat', Orb.T, fmt='%.15E');
np.savetxt('../setup/MC_' + Id_name + '_coef.dat', C.T, fmt='%.15E');

f = open('../setup/MC_' + Id_name + '_config.dat', 'w');
f.write( '%d %d %d %.15f %.15f' % (Npar, Morb, Mdiv, xi, xf) );
f.close();
