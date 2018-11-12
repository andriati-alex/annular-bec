import sys;
import pathlib;
import numpy as np;
import matplotlib as mpl;
import matplotlib.pyplot as plt;

MCTDHBpath = str(pathlib.Path.home()) + '/AndriatiLibrary/annular-bec/pysrc/';
sys.path.insert(0, MCTDHBpath);

import MCTDHBmodule as mc;

suffix_info = '_' + sys.argv[1] + 'time.dat';
nameId = sys.argv[2];

ConfMat = np.loadtxt('../mctdhb_data/' + nameId + '_conf' + suffix_info);



if (ConfMat.ndim > 1):

    NumOfconf = ConfMat.shape[0];

    Morb = int(ConfMat[0][1]);
    Mdx  = int(ConfMat[0][2]);
    xi   = ConfMat[0][3];
    xf   = ConfMat[0][4];
    x    = np.linspace(xi, xf, Mdx + 1);

    for i in range(1, NumOfconf + 1):

        Npar = int(ConfMat[i][0]);

        fname = nameId + 'line-' + str(i) + '_coef' + suffix_info;
        C = np.loadtxt('../mctdhb_data/' + fname, dtype=np.complex128);

        fname = nameId + 'line-' + str(i) + '_orb' + suffix_info;
        orb = np.loadtxt('../mctdhb_data/' + fname, dtype=np.complex128);

        occ = mc.GetOccupation(Npar, Morb, C);
        NO  = mc.GetNatOrb(Npar, Morb, C, orb);
        den = mc.GetGasDensity(occ, NO);
        cor = mc.GetOBcorrelation(occ, NO);

        f = plt.figure(figsize=(8,6));
        ax = plt.gca();
        ax.set_xlim(0.5 * xi, 0.5 * xf);
        ax.set_ylim(0, den.max() * 1.2);
        ax.set_xlabel('Position');
        ax.set_title('Configuration #' + str(i) + ' - Gas Density');
        ax.plot(x, den, x, abs(NO[0])**2);

        f = plt.figure(figsize=(8,6));
        ax = plt.gca();
        Pos = ax.get_position();
        ext = [ xi, xf, xi, xf ];
        im = ax.imshow(1.0 - g, origin='lower', extent=ext, aspect='auto',
                cmap=mpl.cm.magma, interpolation='kaiser');
        cax = f.add_axes( [ Pos.x0 + 1.01 * Pos.width, Pos.y0,
                            0.03 * Pos.width, Pos.height ] )
        f.colorbar(im, ax=ax, cax=cax);
        ax.set_aspect('equal');
        ax.set_xlabel('x');
        ax.set_ylabel('x', rotation=0);
        ax.set_title('Configuration #' + str(i) + ' - One-Body correlation');
        ax.set_aspect('equal');

else :

    i = 1;

    Morb = int(ConfMat[1]);
    Mdx  = int(ConfMat[2]);
    xi   = ConfMat[3];
    xf   = ConfMat[4];
    x    = np.linspace(xi, xf, Mdx + 1);
        
    Npar = int(ConfMat[0]);

    fname = nameId + '_line-' + str(i) + '_coef' + suffix_info;
    C = np.loadtxt('../mctdhb_data/' + fname, dtype=np.complex128);

    fname = nameId + '_line-' + str(i) + '_orb' + suffix_info;
    orb = np.loadtxt('../mctdhb_data/' + fname, dtype=np.complex128);

    occ = mc.GetOccupation(Npar, Morb, C);
    NO  = mc.GetNatOrb(Npar, Morb, C, orb);
    den = mc.GetGasDensity(occ, NO);
    cor = mc.GetOBcorrelation(occ, NO);

    f = plt.figure(figsize=(8,6));
    ax = plt.gca();
    ax.set_xlim(0.5 * xi, 0.5 * xf);
    ax.set_ylim(0, den.max() * 1.2);
    ax.set_xlabel('Position');
    ax.set_title('Gas Density');
    ax.plot(x, den, x, abs(NO[0])**2);

    f = plt.figure(figsize=(8,6));
    ax = plt.gca();
    Pos = ax.get_position();
    ext = [ xi, xf, xi, xf ];
    im = ax.imshow(1.0 - cor, origin='lower', extent=ext,
            aspect='auto', cmap=mpl.cm.magma, interpolation='kaiser');
    ax.set_xlabel('x');
    ax.set_ylabel('x', rotation=0);
    ax.set_title('One-Body correlation');
    cax = f.add_axes( [ Pos.x0 + 0.9 * Pos.width, Pos.y0,
                        0.03 * Pos.width, Pos.height ] )
    f.colorbar(im, ax=ax, cax=cax);
    ax.set_aspect('equal');


plt.show();
