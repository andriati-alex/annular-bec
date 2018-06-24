"""

    SCRIPT MADE TO TEST RUNGE-KUTTA ROUTINE(compare with scipy)

    Before compile the c test program:

    icc -o TestODE -O3 -qopenmp src/tests/TestODE.c -L./lib -I./include -lgp
    
    then run it:

    ./TestODE

    and finally run this script. In red C data points and filled line
    are python data points. Excellent Agreement

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

def f(t, y, params):
    N = y.size;
    lower, mid, upper = params  # unpack parameters
    derivs = np.empty(N, dtype=np.complex128);

    for i in range(1, N - 1):
        derivs[i] = lower[i-1] * y[i-1] + mid[i] * y[i] + upper[i] * y[i+1];

    derivs[0] = mid[0] * y[0] + upper[0] * y[1];
    derivs[N-1] = lower[N-2] * y[N-2] + mid[N-1] * y[N-1];

    return derivs

n = 3;
# Parameters constants that multiply the vector of 
# time dependent functions being solved (RHS of the system)
lower = np.array([  0.8, -1.5 ]);
upper = np.array([ -0.8,  1.5 ]);
mid   = np.array([ -1.4j, -0.5j, 0.2j ]);

# Initial values
y0 = np.array( [ 0.0, 1. / np.sqrt(2.), 1. / np.sqrt(2.) ]);

# Bundle parameters for ODE solver
params = [lower, mid, upper]

# Make time array for solution
t0    = 0.0
tStop = 10.0
tInc  = 0.001
t = np.arange(0., tStop, tInc)

r = ode(f).set_integrator('zvode', atol=1E-10);
r.set_initial_value(y0, t0).set_f_params(params);

# Store all time step solution
S = np.empty( [ t.size, 3 ], dtype=np.complex128);

S[0,:] = y0;
i = 1;
while r.successful() and i < t.size:
    r.integrate(r.t + tInc);
    S[i,:] = r.y;
    i = i + 1;

path = '/home/andriati/AndriatiLibrary/annular-bec/test_out/'
M  = np.loadtxt(path + 'test_ode.dat', dtype=np.complex128);

fig, ax = plt.subplots(1, 3, figsize=(12, 4), sharex=True);

ax[0].plot(t, np.abs(S[:,0]) ** 2);
ax[0].plot(t[::100], np.abs(M[:t.size:100, 0]) ** 2, 'ro', markersize=3);

ax[1].plot(t, np.abs(S[:,1]) ** 2);
ax[1].plot(t[::100], np.abs(M[:t.size:100, 1]) ** 2, 'ro', markersize=3);

ax[2].plot(t, np.abs(S[:,2]) ** 2);
ax[2].plot(t[::100], np.abs(M[:t.size:100, 2]) ** 2, 'ro', markersize=3);

plt.show()
