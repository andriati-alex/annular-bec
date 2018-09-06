"""

    MEMORY USAGE OF MCTDHB METHOD
    *****************************

    Compare the two most known methods to store the coeficients of Fock
    states superposition. The first saves memory but requires much more
    processement(black lines) because requires the mapping of combinato
    rial problem to integers. The  second  goes  exponentially with the
    number of orbitals, with the # of particles as basis,  but requires
    much less processment that the first option(red lines).

"""

from math import factorial as fact;
import numpy as np;
import matplotlib.pyplot as plt;
import matplotlib as mpl;

def f(n, m):
    N = 1;
    for i in range(n + m - 1, n, -1): N = N * i;
    return N / fact(m - 1);

S = np.zeros(1000, dtype=np.int32);
Sexp = np.zeros(1000, dtype=np.int32);

for Npar in range(10, 1001):
    Morb = 2;
    MorbExp = 2;
    while ( f(Npar, Morb) < 1E9 ): Morb = Morb + 1;
    while ( Npar ** (MorbExp - 1) < 1E9 ): MorbExp = MorbExp + 1;
    S[Npar - 1] = Morb - 1;
    Sexp[Npar - 1] = MorbExp - 1;

fig = plt.figure( figsize = (9, 7) );
ax  = plt.axes();

ax.set_xlim(9, 1000);
ax.set_ylim(0, 35);
ax.set_xlabel("# of particles");
ax.set_ylabel("Max. # orbitals");
ax.set_xscale("log");
ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5));

ax.vlines(np.arange(1, 1001, 1)[9:101:5], 0, S[9:101:5]);
ax.vlines(np.arange(1, 1001, 1)[9:101:5], 0, Sexp[9:101:5], 'r');
ax.vlines(np.arange(1, 1001, 1)[101:301:10], 0, S[101:301:10]);
ax.vlines(np.arange(1, 1001, 1)[101:301:10], 0, Sexp[101:301:10], 'r');
ax.vlines(np.arange(1, 1001, 1)[301::20], 0, S[301::20]);
ax.vlines(np.arange(1, 1001, 1)[301::20], 0, Sexp[301::20], 'r');

plt.show();
