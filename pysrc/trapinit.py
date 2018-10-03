
import sys;
import trap;
import numpy as np;

xi   = float(sys.argv[1]); # initial position
xf   = float(sys.argv[2]); # final position
Mdiv = int(sys.argv[3]);   # Number of divisions slices in domain
name = sys.argv[4];        # function to generate Id

x = np.linspace(xi, xf, Mdiv + 1);

f = getattr(trap, name);

params = []; # extra parameters to pass to the function
for i in range(5, len(sys.argv)): params.append( float(sys.argv[i]) );
params = tuple(params);

V = f(x, *params);
np.savetxt('../setup/' + name + '_trap.dat', V.T, fmt='%.15E');
