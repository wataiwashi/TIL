"""
Map
10/12/2024
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import random
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

hide_0 = 1 ### Must be 1
if hide_0 == 1: ### Setting of size, label, etc. 
    if fsize == 0:
        figs = 1.0
        fs1, lw1, ms1 = 1., 1.2, 6.
        tck_s1, alw = 3, 0.625
        ctxt1 = 'paper'
        lpad = [5, 5] # Space between tick label and label
        tpad = [5, 5] # Space between axis and tick label
    else:
        figs = 2.0
        fs1, lw1, ms1 = 1.5, 3., 14.
        tck_s1, alw = 7, 1.25
        ctxt1, ext = 'talk', '_talk' + ext
        lpad = [10, 8]
        tpad = [8, 12]

def sns_set(fs, tck_s, alw, ctxt):
    sns.set(
        context = ctxt,  # Fontsize, linewidth ('paper' or 'talk')
        palette = sns.color_palette("colorblind"),
        font = "serif", 
        font_scale = fs, 
        style = None, 
        rc = {
        'text.usetex': True, 
        'text.latex.preamble' : r'\usepackage{txfonts}\usepackage{bm}',   # LaTeX preamble, txfonts: nature font
        'grid.linestyle': '--', 'grid.linewidth': 0, 
        "xtick.direction":"in", "xtick.major.width":0.8*alw, "xtick.major.size":tck_s, 
        "ytick.direction":"in", "ytick.major.width":0.8*alw, "ytick.major.size":tck_s, 
        "axes.linewidth":alw
        }
    )

def set_ticks_digit(arr_in):
    lis_out = [0]*arr_in.size
    for i in range(arr_in.size):
        j = 0
        lis_out[i] = r'$%d'%(arr_in[i] + 1.0e-6)
        if arr_in[i] + 1.0e-6 < 0.0:
            lis_out[i] = r'$-%d'%(abs(arr_in[i]) + 1.0e-6)
        while (abs(arr_in[i]) + 1.0e-8) % 10.0**j > 1.0e-6:
            if j == 0:
                lis_out[i] = lis_out[i] + r'.'
            lis_out[i] = lis_out[i] + r'%d'%( ( abs(arr_in[i]) % 10.0**(j) )*10.0**(-j + 1) + 1.0e-6 )
            j = j - 1
        lis_out[i] = lis_out[i] + r'$'
    return lis_out
def set_ticks_digit_log(arr_min_max):
    int_log_min_max = [ math.ceil( np.log10(arr_min_max[0]) - 1.0e-6 ), math.floor( np.log10(arr_min_max[-1]) + 1.0e-6 ) ]
    lis_out = []; arr_out = np.empty(0)
    for i in range(int_log_min_max[0], int_log_min_max[1] + 1):
        arr_out = np.append(arr_out, 10**i)
        if i == 0:
            lis_out.append( r'$1$' )  # Set 10^0 as 1
        elif i == 1:
            lis_out.append( r'$10$' ) # Set 10^1 as 10
        else:
            lis_out.append( r'$10^{%d}$'%(i) )
    return arr_out, lis_out

def plot_map():
    fig = plt.figure(figsize = (2.4*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.set_aspect('equal')
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(r'$x_n$', labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(r'$x_{n+1}$', labelpad = lpad[1])
    ax1.tick_params(axis = 'x', pad = tpad[0])
    ax1.tick_params(axis = 'y', pad = tpad[1])
    
    x_map = np.arange(-2.0, 2.0 + 1.0e-6, 1.0e-3)
    y_map = x_map**2

    ax1.plot(x_map, x_map, lw = lw1, ls = 'solid', color = 'C0', alpha = 0.9, zorder = 1, label = r'$x_{n+1}=x_n$')
    ax1.plot(x_map, y_map, lw = lw1, ls = 'solid', color = 'C1', alpha = 0.9, clip_on = False, zorder = 2, label = r'$x_{n+1}=x_n^2$')

    n = 10
    nplt = n*2 - 1
    x0 = np.zeros(nplt)
    x1 = np.zeros(nplt)
    xp0 = np.zeros(nplt)
    xp1 = np.zeros(nplt)
    x0[0] = 0.5
    x1[0] = 1.005
    xp0[0], xp1[0] = -3.0, -3.0
    ax1.scatter(0.0, 0.0, c = 'C2', s = 30, alpha = 0.8, zorder = 3)
    ax1.scatter(1.0, 1.0, c = 'C3', s = 30, alpha = 0.8, zorder = 3)
    for i in range(1, n):
        i1 = i*2 - 1
        i2 = i*2
        xp0[i1], xp1[i1] = x0[i1 - 1]**2, x1[i1 - 1]**2
        xp0[i2], xp1[i2] = x0[i1 - 1]**2, x1[i1 - 1]**2
        x0[i1], x1[i1] = x0[i1 - 1], x1[i1 - 1]
        x0[i2], x1[i2] = xp0[i1], xp1[i1]
    ax1.plot(x0, xp0, lw = lw1, c = 'C2', zorder = 3)
    ax1.plot(x1, xp1, lw = lw1, c = 'C3', zorder = 3)

    ax1.set_xlim(np.amin(x_map), np.amax(x_map)) # Range of axis
    ax1.set_ylim(np.amin(x_map)*0.5, np.amax(y_map))
    ax1.set_xticks(ax1.get_xticks(), set_ticks_digit(ax1.get_xticks()))
    ax1.set_yticks(ax1.get_yticks(), set_ticks_digit(ax1.get_yticks()))

    ### Legend
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ## Save file
    fig.savefig(plotdir + "map" + ext, bbox_inches = "tight")

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    sns_set(fs1, tck_s1, alw, ctxt1)

    plot_map()

    if fshow == 1:
        plt.show()
    print("end main")
