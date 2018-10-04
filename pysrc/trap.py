
import numpy as np;

def zero(x): return np.zeros(x.size);





def harmonic(x, omega): return (omega * omega / 2) * x * x;





def deltaBarrier(x, height):
    V = np.zeros(x.size, dtype=np.float128);
    V[int((x.size - 1) / 2)] = height / (x[1] - x[0]);
    return V;





def doubleDeltaBarrier(x, h1, h2):
    V = np.zeros(x.size, dtype=np.float128);
    V[x.size / 2] = h1 / (x[1] - x[0]);
    V[x.size - 1] = h2 / (x[1] - x[0]);
    V[0] = h2 / (x[1] - x[0]);
    return V;





def DoubleWell(x, omega, a):
    return (omega**2 / 2) * (x**2 - a**2)**2;





def narrowDoubleWell(x, omega, a):
    return (omega**2 / 2) * (x**2 - a**2)**2 / a**5;





def tightDoubleWell(x, omega, a):
    return (omega**2 / 2) * (x**2 - a**2)**2 / a**4;





def deltaDoubleWell(x, omega, h):
    V = (omega * omega / 2) * x * x;
    V[x.size / 2] = h / (x[1] - x[0]);
    return V;
