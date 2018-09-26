import sys;
import numpy as np;
from scipy.integrate import simps;

"""

    GENERATE INITIAL CONDITION DATA
    *******************************

    Given a keyword generate one of some pre-defined initial condition
    data, like solitons, or others yet to be implemented. This initial
    condition can be take to be time evolved or as initial guess  to a
    method to obtain stationary solution.

    CALL
    ****

    python generate_init.py x1 x2 M Id extra_params

    where [x1, x2] is the position domain and M the number of
    discretized intervals. The number of points therefore  is
    M + 1 in an vector holding each position value. Id is the
    function Identification to generate data and extra_params
    to be able to evalute it.

"""

lf = np.float128;
lc = np.complex256;

def FourierLocModes(x, a, c, n):
    S = BrightSoliton(x, a, c);
    N = np.sqrt( simps(abs(S)**2, dx=(x[1]-x[0])) );
    # generate random numbers in the range [-1, 1]
    noise = (np.random.random(int(n)) - 0.5) / 0.5;
    # Localized by a Gaussian like-shape
    sig = 0.2 * (x[-1] - x[0]);
    mid = (x[-1] + x[0]) / 2;
    # Fourier modes with some noise
    k = (np.arange(1, int(n) + 1) + noise) * 2 * np.pi / sig;
    weights = 1.0 / np.arange(1, n + 1);
    for i in range(int(n)):
        S = S + weights[i] * np.exp(1.0j * k[i] * x - ((x - mid) / sig) ** 2);
    return S * N / np.sqrt( simps(abs(S)**2, dx=(x[1]-x[0])) );

def BrightSoliton(x, a, c):
    numerator = a * np.exp(0.5j * c * x * np.sqrt(2), dtype=lc);
    denominator = np.cosh(a * x / np.sqrt(2), dtype=lf);
    return numerator / denominator;

def NBrightSoliton(x, a, c):
    numerator = a * np.exp(0.5j * c * x * np.sqrt(2), dtype=lc);
    denominator = np.cosh(a * x / np.sqrt(2), dtype=lf);
    return numerator / denominator / np.sqrt(2 * np.sqrt(2) * a);

def DarkSolitonModes(x, a, b, c, n):
    S = np.sqrt(a + b * (np.tanh(x / c) ** 2));
    N = np.sqrt( simps(abs(S)**2, dx=(x[1]-x[0])) );
    # periodic frequency modes
    k = np.arange(1, n + 1) * 2 * np.pi / (x[-1] - x[0]);
    t = (np.random.random(int(n)) - 0.5) / 3;
    weights = t * b / ( np.arange(1, int(n) + 1) );
    for i in range(int(n)):
        S = S + weights[i] * np.exp(1.0j * k[i] * x);
    return S * N / np.sqrt( simps(abs(S)**2, dx=(x[1]-x[0])) );

# Domain discretization parameters
x1 = lf(sys.argv[1]);
x2 = lf(sys.argv[2]);
M  = int(sys.argv[3]);
Id = int(sys.argv[4]);

Params = [];
for i in range(5, len(sys.argv)): Params.append(lf(sys.argv[i]));
Params = tuple(Params);

x  = np.linspace(x1, x2, M + 1, dtype=lf);
dx = (x2 - x1) / M;

out = np.empty(x.size, dtype=lc);

if (Id == 1):
    out = NBrightSoliton(x, *Params);
    Id_name = 'BrightSoliton';
elif (Id == 2):
    out = FourierLocModes(x, *Params)
    Id_name = 'NoiseBright';
elif (Id == 3):
    out = DarkSolitonModes(x, *Params)
    Id_name = 'NoiseDark';

np.savetxt('setup/' + Id_name + '_init.dat', out, fmt='%.15E');

f = open('setup/' + Id_name + '_domain.dat', "w");

f.write("%.2f %.2f %d" % (x1, x2, M));
f.close();
