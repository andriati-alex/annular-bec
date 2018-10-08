
import sys;
import trap;
import numpy as np;
from pathlib import Path;

"""
===========================================================================

SCRIPT TO GENERATE (TRAP)POTENTIAL IN DISCRETIZED POSITIONS
-----------------------------------------------------------

    Calling:
    -------
    python trapinit.py xi xf Mdiv Potential_name par1 ... parN out_filename

    Description:
    -----------
    Call an implemented routine in trap.py module that return a numpy array
    with values of the potential in each discretized potition

===========================================================================
"""

xi   = float(sys.argv[1]); # initial position
xf   = float(sys.argv[2]); # final position
Mdiv = int(sys.argv[3]);   # Number of divisions slices in domain
name = sys.argv[4];        # function to generate Id

x = np.linspace(xi, xf, Mdiv + 1);

f = getattr(trap, name);

params = []; # extra parameters to generate potential values
for i in range(5, len(sys.argv) - 1): params.append( float(sys.argv[i]) );
params = tuple(params);

outname = sys.argv[len(sys.argv) - 1];

home = str(Path.home());

outname = home + '/AndriatiLibrary/annular-bec/setup/' + outname

V = f(x, *params);
np.savetxt(outname + '_trap.dat', V.T, fmt='%.15E');
