import sys;
import numpy as np;
import matplotlib as mpl;
from pathlib import Path;

"""

    SHOW RESULT OF TIME EVOLUTION IN ANIMATION
    ******************************************

    Given a file ID after calling time_evolution executable show the
    density |f(x,t)| ^ 2 where each fstep is a frame yeilding a movie
    of time evolution.

    CALL
    ****

    $ python watch_evolution.py file_id frame_fstep

    file_id    - A valid name after running time_evolution
    frame_fstep - Jump some time steps of solution to show a shorter movie

"""

# IMPORT PLOT MODULES #
import matplotlib.pyplot as plt;
from matplotlib.animation import FuncAnimation;

folder = str(Path.home()) + '/AndriatiLibrary/annular-bec/output/';

fname = folder + sys.argv[1] + '_line-1_orb_realtime.dat';
fstep = int(sys.argv[2]); # how much time-steps a frame jumps

S = np.loadtxt(fname, dtype=np.complex128);
Smod2 = np.absolute(S)[::fstep,:] ** 2;

fname = folder + sys.argv[1] + '_conf_realtime.dat'
domain = np.loadtxt(fname, dtype=np.float64);

y2 = Smod2.max() + 0.1 * (Smod2.max() - Smod2.min());
y1 = Smod2.min() - 0.1 * (Smod2.max() - Smod2.min());

x1 = domain[1];
x2 = domain[2];
grid = int(domain[0]);
x = np.linspace(x1, x2, grid);





fig   = plt.figure();
ax    = plt.axes();
line, = ax.plot([], [], '-');

def init():
    line.set_data([], []);
    ax.set_xlim(x1,x2);
    ax.set_ylim(y1,y2);
    return line,

def update(frame):
    line.set_data(x, Smod2[frame]);
    return line,

frames = np.arange(Smod2.shape[0], dtype=np.int32);

anim = FuncAnimation(fig, update, init_func=init, frames=frames,
       interval=20, blit=True);

plt.show();
