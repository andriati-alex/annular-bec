import sys;
import numpy as np;
import matplotlib as mpl;

"""

    SHOW RESULT OF TIME EVOLUTION IN ANIMATION
    ******************************************

    Given a file ID after calling time_evolution executable show the
    density |f(x,t)| ^ 2 where each step is a frame yeilding a movie
    of time evolution.

    CALL
    ****

    $ python watch_evolution.py file_id frame_step

    file_id    - A valid name after running time_evolution
    frame_step - Jump some time steps of solution to show a shorter movie

"""

# SETUP PLOT PARAMS   #
mpl.use('Qt5Agg');
# IMPORT PLOT MODULES #
import matplotlib.pyplot as plt;
from matplotlib import animation;

fname = sys.argv[1];
step  = int(sys.argv[2]);

S = np.loadtxt('../gp_data/' + fname + '_time.dat', dtype=np.complex128);
Smod2 = np.absolute(S)[::step,:] ** 2;

domain = np.loadtxt('../gp_data/' + fname + '_domain.dat', dtype=np.float64);

y2 = Smod2.max() + 0.1 * (Smod2.max() - Smod2.min());
y1 = Smod2.min() - 0.1 * (Smod2.max() - Smod2.min());

x1 = domain[0];
x2 = domain[1];
x  = np.linspace(x1, x2, S.shape[1]);

fig   = plt.figure(figsize=(10, 8));
ax    = plt.axes(xlim=(x1, x2), ylim=(y1, y2));
line, = ax.plot([], [], lw=2);

def init():
    line.set_data([], []);
    return line,

def anim_frame(i):
    line.set_data(x, Smod2[i,:]);
    return line,

anim = animation.FuncAnimation(fig, anim_frame, init_func=init,
       frames=Smod2.shape[0], interval=20, blit=True);

plt.show();
