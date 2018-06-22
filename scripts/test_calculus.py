"""

    SCRIPT MADE TO TEST INTEGRATION AND DERIVATIVES

    f(x) = cos(2 * PI * x / 4) + I * sin(PI * x / 4 - PI / 2);

    The plot using circle marks are obtained using finite differences
    to compute the derivatives.

    Before compile the c test program:

    icc -o TestCalculus -O3 -qopenmp -mkl src/tests/TestCalculus.c -L./lib
    -I./include -lgp
    
    then run it:

    ./TestCalculus

"""

import numpy as np
import matplotlib.pyplot as plt

path = '/home/andriati/AndriatiLibrary/linear-algebra/test_out/'
M = np.loadtxt(path + 'calculus_out.dat', dtype=np.complex128);

x = np.linspace(-8, 8, M.shape[1]);
f = M[0,:];
dfdx = M[1,:];
dfdxdiff = M[2,:];
ints = M[3,:];

fig, ax = plt.subplots(1, 3, figsize=(14, 4), sharex=True);

ax[0].plot(x, f.real, 'b');
ax[0].plot(x, f.imag, 'r');
ax[0].set_xlabel("x");
ax[0].set_title("f");

ax[1].plot(x, dfdx.real, 'b');
ax[1].plot(x[::200], dfdxdiff[::200].real, 'ko', markersize=2.0);
ax[1].plot(x, dfdx.imag, 'r');
ax[1].plot(x[::400], dfdxdiff[::400].imag, 'ko', markersize=2.0);
ax[1].set_xlabel("x");
ax[1].set_title("dfdx");

ax[2].plot(x, ints.real, 'b');
ax[2].plot(x, ints.imag, 'r');
ax[2].set_xlabel("x");
ax[2].set_title("integral of f");

plt.show()
