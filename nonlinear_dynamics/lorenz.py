"""
Lorenz
08/04/2024
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

fshow = 1 # Display popup 0: no, 1: yes
### Output setting
fsize = 0   # 0: small size for paper, 1: large size for talk
ext = 'pdf' # File extension (pdf, svg, png, etc.)
plotdir = 'plot' # Plot file save directry

datadir = 'data'
plotdir, ext = plotdir + '/', '.' + ext
dfile   = datadir + 'data1.dat'
os.makedirs(plotdir, exist_ok = True) # Make plotdir
if os.path.exists(plotdir) == False:
    sys.exit("plotdir error: %s"%plotdir)

### Lorenz parameter
sigma = 10.0
b = 8.0/3.0
r = 28.0
# r = 0.8
dt = 1.0e-4
tmax, dt_plt = 50.0, 0.5e-2

def lorenz(x0, y0, z0, tmax, dt_plt):
    ntmax = int((tmax + 1.0e-8)/dt) + 1
    n_p = int((tmax + 1.0e-8)/dt_plt) + 1
    t_plt = np.zeros(n_p)
    xyz_plt = np.zeros((3, n_p))
    nt_plt = 0
    xyz = np.array([x0, y0, z0])
    dot_xyz = np.zeros(3)
    for nt in range(ntmax):
        t = nt*dt
        if nt == nt_plt:
            i_p = int((t + 1.0e-8)/dt_plt)
            t_plt[i_p] = t
            xyz_plt[:, i_p] = xyz[:]
            nt_plt = int((t + dt_plt + 1.0e-8)/dt)
        dot_xyz[0] = -sigma*xyz[0] + sigma*xyz[1]
        dot_xyz[1] = -xyz[0]*xyz[2] + r*xyz[0] - xyz[1]
        dot_xyz[2] = xyz[0]*xyz[1] - b*xyz[2]
        xyz[:] = xyz[:] + dot_xyz[:]*dt
    return t_plt, xyz_plt

t_plt1, xyz_plt1 = lorenz(0.0, 0.1, 0.0, tmax, dt_plt)
t_plt2, xyz_plt2 = lorenz(0.0, 0.10001, 0.0, tmax, dt_plt)

hide_0 = 1 ### Must be 1
if hide_0 == 1: ### Setting of size, label, etc. 
    lx, ly, lz = r'$X$', r'$Y$', r'$Z$' # Use TeX character with r''
    if fsize == 0:
        figs = 1.
        fs1, lw1, ms1 = 1., 1.2, 6.
        tck_s1, alw = 3, 0.625
        ctxt1 = 'paper'
        lpad = [5, 5] # Space between tick label and label
        tpad = [3, 5] # Space between axis and tick label
    else:
        figs = 2.
        fs1, lw1, ms1 = 1.5, 3., 14.
        tck_s1, alw = 7, 1.25
        ctxt1, ext = 'talk', '_talk' + ext
        lpad = [10, 8]
        tpad = [8, 12]

def sns_set(fs, tck_s, alw, ctxt):
    sns.set(
        context = ctxt,  # Fontsize, linewidth ('paper' or 'talk')
        palette = sns.color_palette("colorblind"),
        font = "serif",       # Font
        font_scale = fs,  # Font scale (Changing this to further tweak preset determined by context)
        style = None, 
        # style = "whitegrid",   # White background with grid
        rc = {
        # 'text.usetex': True, 
        # 'text.latex.preamble' : r'\usepackage{txfonts}',   # LaTeX preamble, txfonts: nature font
        'grid.linestyle': '--', 'grid.linewidth': 0, 
        "xtick.direction":"in", "xtick.major.width":0.8*alw, "xtick.major.size":tck_s, 
        "ytick.direction":"in", "ytick.major.width":0.8*alw, "ytick.major.size":tck_s, 
        "axes.linewidth":alw
        }
    )

def plot_y():
    fig = plt.figure(figsize = (4.8*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lx, labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(ly, labelpad = lpad[1])
    xm = [0.0, tmax]
    ax1.set_xlim(xm[0], xm[1]) # Range of axis
    ax1.tick_params(axis = 'x', pad = tpad[0])
    ax1.tick_params(axis = 'y', pad = tpad[1])

    ax1.plot(t_plt1, xyz_plt1[1], lw = lw1, ls = 'solid', color = 'C0', alpha = 1.0, clip_on = False, zorder = 8, label = r'$Y(t=0)=%8.5f$'%xyz_plt1[1, 0])

    ax1.plot(t_plt2, xyz_plt2[1], lw = lw1, ls = 'solid', color = 'C1', alpha = 0.8, clip_on = False, zorder = 8, label = r'$Y(t=0)=%8.5f$'%xyz_plt2[1, 0])

    ### Legend
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ## Save file
    fig.savefig(plotdir + "lorenz_y" + ext, bbox_inches = "tight")

def plot_xyz():
    fig = plt.figure(figsize = (3.6*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111, projection = '3d')
    ax1.set_xlabel(lx) # Axis label
    ax1.set_ylabel(ly)
    ax1.set_zlabel(lz)

    ax1.plot(xyz_plt1[0], xyz_plt1[1], xyz_plt1[2], lw = lw1*0.25, ls = 'solid', color = 'C0', alpha = 1.0, clip_on = False, zorder = 8, label = r'$y_0=%10.5f$'%xyz_plt1[1, 0])
    # ax1.plot(xyz_plt2[0], xyz_plt2[1], xyz_plt2[2], lw = lw1, ls = 'solid', color = 'C1', alpha = 0.8, clip_on = False, zorder = 8, label = r'$y_0=%10.5f$'%xyz_plt2[1, 0])

    ### Legend
    # h1, l1 = ax1.get_legend_handles_labels()
    # ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ### save
    fig.savefig(plotdir + 'lorenz.jpg', dpi = 600, bbox_inches = 'tight', pad_inches = 0.4)

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("plotdir: ", plotdir, ", ctxt: ", ctxt1)
    print('ntmax, n_p: ', int((tmax + 1.0e-8)/dt), int((tmax + 1.0e-8)/dt_plt))
    sns_set(fs1, tck_s1, alw, ctxt1)

    plot_y()
    plot_xyz()

    if fshow == 1:
        plt.show()
    print("end main")
