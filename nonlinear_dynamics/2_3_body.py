"""
2 body problem and 3 body problem
09/04/2024
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

fshow = 1 # Display popup 0: no, 1: yes
### Output setting
fsize = 0   # 0: small size for paper, 1: large size for talk
ext = 'svg' # File extension (pdf, svg, png, etc.)
plotdir = 'plot' # Plot file save directry

datadir = 'data'
plotdir, ext = plotdir + '/', '.' + ext
dfile   = datadir + 'data1.dat'
os.makedirs(plotdir, exist_ok = True) # Make plotdir
if os.path.exists(plotdir) == False:
    sys.exit("plotdir error: %s"%plotdir)

### parameter
G = 10.0
m = np.array([1.0, 1.0, 0.1])
dt = 1.0e-4
tmax, dt_plt = np.array([150.0, 8.0]), np.array([2.0e0, 5.0e-2])
### initial condition
x1_0 = np.array([-1.0,  0.0])
x2_0 = np.array([ 1.0,  0.0])
x3_0 = np.array([ 0.0,  3.0])
v1_0 = np.array([ 0.0, -1.0])
v2_0 = np.array([ 0.0,  1.0])
v3_0 = np.array([-0.5,  0.5])

def two_body(x1_0, x2_0, v1_0, v2_0, tmax, dt_plt):
    ntmax = int((tmax + 1.0e-8)/dt) + 1
    n_p = int((tmax + 1.0e-8)/dt_plt) + 1
    t_plt = np.zeros(n_p)
    x1_plt, x2_plt = np.zeros((2, n_p)), np.zeros((2, n_p))
    nt_plt = 0
    x1, x2 = x1_0, x2_0
    v1, v2 = v1_0, v2_0
    for nt in range(ntmax):
        t = nt*dt
        if nt == nt_plt:
            i_p = int((t + 1.0e-8)/dt_plt)
            t_plt[i_p] = t
            x1_plt[:, i_p], x2_plt[:, i_p] = x1[:], x2[:]
            nt_plt = int((t + dt_plt + 1.0e-8)/dt)
        r12 = np.sqrt( (x2[0] - x1[0])**2 + (x2[1] - x1[1])**2 )
        v1 = v1 - G*m[0]*m[1]*(x1 - x2)/r12**3*dt/m[0]
        v2 = v2 - G*m[0]*m[1]*(x2 - x1)/r12**3*dt/m[1]
        x1, x2 = x1 + v1*dt, x2 + v2*dt
    return t_plt, x1_plt, x2_plt

def three_body(x1_0, x2_0, x3_0, v1_0, v2_0, v3_0, tmax, dt_plt):
    ntmax = int((tmax + 1.0e-8)/dt) + 1
    n_p = int((tmax + 1.0e-8)/dt_plt) + 1
    t_plt = np.zeros(n_p)
    x1_plt, x2_plt, x3_plt = np.zeros((2, n_p)), np.zeros((2, n_p)), np.zeros((2, n_p))
    nt_plt = 0
    x1, x2, x3 = x1_0, x2_0, x3_0
    v1, v2, v3 = v1_0, v2_0, v3_0
    for nt in range(ntmax):
        t = nt*dt
        if nt == nt_plt:
            i_p = int((t + 1.0e-8)/dt_plt)
            t_plt[i_p] = t
            x1_plt[:, i_p], x2_plt[:, i_p], x3_plt[:, i_p] = x1[:], x2[:], x3[:]
            nt_plt = int((t + dt_plt + 1.0e-8)/dt)
        r12 = np.sqrt( (x2[0] - x1[0])**2 + (x2[1] - x1[1])**2 )
        r23 = np.sqrt( (x2[0] - x3[0])**2 + (x2[1] - x3[1])**2 )
        r31 = np.sqrt( (x3[0] - x1[0])**2 + (x3[1] - x1[1])**2 )
        v1 = v1 - G*( m[0]*m[1]*(x1 - x2)/r12**3 + m[0]*m[2]*(x1 - x3)/r31**3 )*dt/m[0]
        v2 = v2 - G*( m[1]*m[0]*(x2 - x1)/r12**3 + m[1]*m[2]*(x2 - x3)/r23**3 )*dt/m[1]
        v3 = v3 - G*( m[2]*m[0]*(x3 - x1)/r31**3 + m[2]*m[1]*(x3 - x2)/r23**3 )*dt/m[2]
        x1, x2, x3 = x1 + v1*dt, x2 + v2*dt, x3 + v3*dt
    return t_plt, x1_plt, x2_plt, x3_plt

hide_0 = 1 ### Must be 1
if hide_0 == 1: ### Setting of size, label, etc. 
    lt, lx, ly = r'$t$', r'$x$', r'$y$' # Use TeX character with r''
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

def plot_two_body():
    print('ntmax, n_p: ', int((tmax[0] + 1.0e-8)/dt), int((tmax[0] + 1.0e-8)/dt_plt[0]))
    t_plt1, x1_plt, x2_plt = two_body(x1_0, x2_0, v1_0, v2_0, tmax[0], dt_plt[0])

    fig = plt.figure(figsize = (4.8*figs, 3.6*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw); ax1.spines["left"].set_linewidth(alw); ax1.spines["bottom"].set_linewidth(alw); ax1.spines["right"].set_linewidth(alw)
    ax1.set_aspect ('equal')

    ax1.set_xlabel(lx, labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(ly, labelpad = lpad[1])
    ax1.tick_params(axis = 'both', pad = tpad[0])

    ax1.scatter(x1_plt[0, 0], x1_plt[1, 0], s = 50, marker = '*', c = 'red', zorder = 10)
    ax1.scatter(x2_plt[0, 0], x2_plt[1, 0], s = 50, marker = '*', c = 'red', zorder = 10)
    ax1.scatter(x1_plt[0], x1_plt[1], marker = 'o', color = 'C0', alpha = 1.0, clip_on = False, zorder = 8, label = r'1')
    ax1.scatter(x2_plt[0], x2_plt[1], marker = 'o', color = 'C1', alpha = 0.8, clip_on = False, zorder = 8, label = r'2')

    ### Legend
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ## Save file
    fig.savefig(plotdir + "2_body" + ext, bbox_inches = "tight")

def plot_three_body():
    print('ntmax, n_p: ', int((tmax[1] + 1.0e-8)/dt), int((tmax[1] + 1.0e-8)/dt_plt[1]))
    t_plt1, x1_plt, x2_plt, x3_plt = three_body(x1_0, x2_0, x3_0, v1_0, v2_0, v3_0, tmax[1], dt_plt[1])

    fig = plt.figure(figsize = (4.8*figs, 3.6*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw); ax1.spines["left"].set_linewidth(alw); ax1.spines["bottom"].set_linewidth(alw); ax1.spines["right"].set_linewidth(alw)
    ax1.set_aspect ('equal')

    ax1.set_xlabel(lx, labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(ly, labelpad = lpad[1])
    ax1.tick_params(axis = 'both', pad = tpad[0])

    ax1.scatter(x1_plt[0, 0], x1_plt[1, 0], s = 50, marker = '*', c = 'red', zorder = 10)
    ax1.scatter(x2_plt[0, 0], x2_plt[1, 0], s = 50, marker = '*', c = 'red', zorder = 10)
    ax1.scatter(x3_plt[0, 0], x3_plt[1, 0], s = 50, marker = '*', c = 'red', zorder = 10)
    ax1.scatter(x1_plt[0], x1_plt[1], marker = 'o', color = 'C0', alpha = 1.0, clip_on = False, zorder = 8, label = r'1')
    ax1.scatter(x2_plt[0], x2_plt[1], marker = 'o', color = 'C1', alpha = 0.8, clip_on = False, zorder = 8, label = r'2')
    ax1.scatter(x3_plt[0], x3_plt[1], marker = 'o', color = 'C2', alpha = 0.6, clip_on = False, zorder = 8, label = r'3')

    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    fig.savefig(plotdir + "3_body" + ext, bbox_inches = "tight")

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)

    plot_two_body()
    plot_three_body()

    if fshow == 1:
        plt.show()
    print("end main")
