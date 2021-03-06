
import sys;
import numpy as np;
from scipy.integrate import simps;
from pathlib import Path;





"""
=============================================================================


    GENERATE INITIAL CONDITION DATA
    -------------------------------

    Given a keyword generate one of some pre-defined initial condition
    data, like solitons, or others yet to be implemented. This initial
    condition can be take to be time evolved or as initial guess  to a
    method to obtain stationary solution for Gross-Pitaevskii equation



    CALL IN COMMAND LINE
    --------------------

    python generate_init.py x1 x2 M Id extra_params

    where [x1, x2] is the position domain and M the number of
    discretized intervals. The number of points therefore  is
    M + 1 in an vector holding each position value. Id is the
    function Identification to generate data and extra_params
    to be able to evalute it.


=============================================================================
"""





lf = np.float128;
lc = np.complex256;





def HarmonicTrap(x, a, n):
    n = int(n);
    # generate random numbers in the range [-0.5, 0.5]
    noise = (np.random.random(int(n)) - 0.5);
    # Localized by a Gaussian like-shape
    sig = 0.14 * (x[-1] - x[0]);
    mid = (x[-1] + x[0]) / 2;
    # Fourier modes with some noise
    S = BrightSoliton(x, a, 1.0 / sig);
    k = (np.arange(0, int(n)) + noise) * 2 * np.pi / sig;
    weights = a / np.arange(1, n + 1)**1.5;
    for i in range(int(n)):
        plus  = +1.0j * k[i] * x - ((x - mid) / sig)**2;
        minus = -1.0j * k[i] * x - ((x - mid) / sig)**2;
        S = S + noise[i] * weights[i] * np.exp(plus);
        S = S + noise[n-1-i] * weights[i] * np.exp(minus);
    return S / np.sqrt( simps( abs(S)**2, dx = x[1] - x[0] ) );





def BrightSoliton(x, a, c):
    numerator = a * np.exp(0.5j * c * x * np.sqrt(2), dtype=lc);
    denominator = np.cosh(a * x / np.sqrt(2), dtype=lf);
    return numerator / denominator;





def NormBrightSoliton(x, a, c):
    BS = BrightSoliton(x, a, c);
    return BS / np.sqrt( simps( abs(BS)**2, dx = x[1] - x[0] ) );





def Ring(x, a, b, n):
    n = int(n);
    k = np.arange(1, n + 1) * 2 * np.pi / (x[-1] - x[0]);
    t = (np.random.random(int(n)) - 0.5) / 0.5;
    c = (x[-1] - x[0]) / 10;
    S = np.sqrt(a + b * (np.tanh(x / c) ** 2)) * np.exp(2 * np.pi * t[0]);
    shape = np.sqrt(a + b * (np.tanh(x / c) ** 2));
    weights = 0.5 * b / ( np.arange(1, int(n) + 1) );
    for i in range(int(n)):
        plus  =  1.0j * (k[i] * x);
        minus = -1.0j * (k[i] * x);
        S = S + t[i] * shape * weights[i] * np.exp(plus);
        S = S + t[n-1-i] * shape * weights[i] * np.exp(minus);
    return S / np.sqrt( simps(abs(S)**2, dx = x[1] - x[0] ) );





def RingBarrier(x, a, omega, n):
    n = int(n);
    c = (x[-1] - x[0]) / 10;
    freq = 2 * np.pi / (x[-1] - x[0]);
    t = (np.random.random(int(n)) - 0.5) / 0.5;
    S = (1 - np.exp( - abs(x / c) + 1.0j * x * freq) / a);
    for m in range(1, n):
        k = m + 1;
        S = S + np.exp(1.0j*freq*omega*x)*np.sin(m*freq*x)*t[m]/(k*k);
    return S / np.sqrt( simps(abs(S)**2, dx = x[1] - x[0] ) );





#              Domain discretization parameters               #
#              --------------------------------               #





x1 = lf(sys.argv[1]);
x2 = lf(sys.argv[2]);
M = int(sys.argv[3]);
fname = sys.argv[4];

Params = [];
for i in range(5, len(sys.argv)): Params.append(lf(sys.argv[i]));
Params = tuple(Params);

x  = np.linspace(x1, x2, M + 1, dtype=lf);
dx = (x2 - x1) / M;





if   (fname == 'BrightSoliton') :     out = BrightSoliton(x, *Params);
elif (fname == 'NormBrightSoliton') : out = NormBrightSoliton(x, *Params);
elif (fname == 'HarmonicTrap') :      out = HarmonicTrap(x, *Params)
elif (fname == 'Ring') :              out = Ring(x, *Params)
elif (fname == 'RingBarrier') :       out = RingBarrier(x, *Params)
else : raise IOError('\n\nInitial function name not implemented.\n\n');





folder = str(Path.home()) + '/AndriatiLibrary/annular-bec/input/';

np.savetxt(folder + fname + '_init.dat', out.T, fmt='%.15E');

f = open(folder + fname + '_domain.dat', "w");

f.write("%.15f %.15f %d" % (x1, x2, M));

f.close();
