import sys;
import numpy as np;
from scipy.integrate import simps;
from math import pi;
from math import sqrt;
from math import factorial as fac;

"""

    GENERATE INITIAL CONDITION DATA TO MCTDHB
    *****************************************

    Generate 3 out files in setup/ folder containing:

    1) _config.dat

        (# of particles) (# of orbitals) (# of divisions in x)

    2) _orb.dat

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

    3) _coef.dat

        A column vector of length:

             (Npar + Morb - 1)!
            --------------------
            (Npar)!  (Morb - 1)!

"""

lf = np.float128;
lc = np.complex256;

def NC(Npar, Morb):
    n = 1;
    if (Morb > Npar):
        for i in range(Npar + Morb - 1, Morb - 1, -1): n = n * i;
        return int(n / fac(Npar));
    else:
        for i in range(Npar + Morb - 1, Npar, -1): n = n * i;
        return int(n / fac(Morb - 1));

def IndexToFock(k, N, M, v):
    x = 0;
    m = M - 1;
    for i in range(M): v[i] = 0;
    while (k > 0):
        while (k - NC(N, m) < 0): m = m - 1;
        x = k - NC(N, m);
        while (x >= 0):
            v[m] = v[m] + 1;
            N = N - 1;
            k = x;
            x = x - NC(N, m);
    for i in range(N, 0, -1): v[0] = v[0] + 1;

def renormalize(f, dx):
    renorm = sqrt(1.0 / simps(abs(f)**2, dx=dx));
    return f * renorm;

def BrightSoliton(x, a, c):
    numerator = a * np.exp(0.5j * c * x * np.sqrt(2), dtype=lc);
    denominator = np.cosh(a * x / np.sqrt(2), dtype=lf);
    return numerator / denominator / np.sqrt(2 * np.sqrt(2) * a);

def AngularMom(Morb, x, S):
    S[0,:] = 1.0 / sqrt(2 * pi);
    for i in range(2, Morb, 2):
        S[i - 1,:] = np.exp(- 1.0j * (i / 2) * x, dtype=lc) / sqrt(2 * pi);
        S[i, :] = np.exp(1.0j * (i / 2) * x, dtype=lc) / sqrt(2 * pi);

def NoiseAngular(Morb, x, S, k):
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

def Coef(Npar, Morb, C):
    v = np.empty(Morb, dtype=np.int32);
    phase = np.exp(2 * pi * 1.0j * (np.random.random(NC(Npar, Morb)) - 0.5));
    for l in range(NC(Npar, Morb)):
        IndexToFock(l, Npar, Morb, v);
        prod = 1.0;
        for j in range(Morb):
            if (j % 2) == 0: k = j - 1;
            else           : k = j;
            prod = prod * np.exp( - float(v[j] * k * k) / Npar, dtype=lf);
        C[l] = phase[l] * prod;
    Norm = 0;
    for l in range(NC(Npar, Morb)): Norm = Norm + abs(C[l]) * abs(C[l]);
    for l in range(NC(Npar, Morb)): C[l] = C[l] / np.sqrt(Norm);

Npar = int(sys.argv[1]); # Number of Particles
Morb = int(sys.argv[2]); # Number of orbitals
Mdiv = int(sys.argv[3]); # Number of divisions between -pi to pi

x  = np.linspace(-8, 8, Mdiv + 1, dtype=lf);
dx = 2 * 8.0 / Mdiv;

Orb = np.zeros([Morb, x.size], dtype=lc); # orbitals
C = np.zeros(NC(Npar, Morb), dtype=lc);   # coeficients

C[0] = 1;
Orb[0] = BrightSoliton(x, 2.0, 3.0);
# Coef(Npar, Morb, C);

Id_name = 'BrightSoliton-' + str(Npar) + '-' + str(Morb);

np.savetxt('setup/MC_' + Id_name + '_orb.dat', Orb.T, fmt='%.15E');
np.savetxt('setup/MC_' + Id_name + '_coef.dat', C.T, fmt='%.15E');

f = open('setup/MC_' + Id_name + '_config.dat', "w");
f.write("%d %d %d %.15f %.15f" % (Npar, Morb, Mdiv, -8, 8));
f.close();
