"""
Log plot etc. LaTeX font
12/02/2024
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import os
import sys

fshow = 0   # Display popup 0: no, 1: yes
### Output setting
fsize = 0   # 0: small size for paper, 1: large size for talk
ext = 'pdf' # File extension (pdf, svg, png, etc.)
plotdir = 'plot' # Plot file save directry

datadir = 'data'
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
dfile   = datadir + 'data1.dat'

if os.path.exists(datadir) == False:  # Cancel if datadir or plotdir does not exist
    sys.exit("datadir error: %s"%datadir)
os.makedirs(plotdir, exist_ok = True) # Make plotdir
if os.path.exists(plotdir) == False:
    sys.exit("plotdir error: %s"%plotdir)

### Load text file
x1 = np.loadtxt(dfile, usecols = 0, dtype = 'float32') # usecols: column number, dtype (float32: single precision, float64: double precision, int64: integer, etc.)
y3 = np.loadtxt(dfile, usecols = 3, dtype = 'float32')

### Setting of size, label, etc.
lx, ly = r'Distance, $x$', r'$y$' # Use TeX character with r''
color0 = 'C0'
color2 = 'C2'
if fsize == 0:
    figs = 1.
    fs1, lw1, ms1 = 1., 1., 6.
    tck_s1, alw = 3, 0.625
    ctxt1 = 'paper'
    lpad = [5, 5] # Space between tick label and label
    tpad = 8      # Space between axis and tick label
else:
    figs = 2.
    fs1, lw1, ms1 = 1.5, 3., 14.
    tck_s1, alw = 7, 1.25
    ctxt1, ext = 'talk', '_talk' + ext
    lpad = [10, 8]
    tpad = 15

def sns_set(fs, tck_s, alw, ctxt):
    sns.set(
        context = ctxt,      # Fontsize, linewidth ('paper' or 'talk')
        palette = sns.color_palette("colorblind"),
        font = "serif",      # Font
        # font = "sans-serif",
        font_scale = fs,     # Font scale (Changing this to further tweak preset determined by context)
        style = None, 
        # style = "whitegrid", # White background with grid
        rc = {
        'text.usetex': True, 
        'text.latex.preamble' : r'\usepackage{txfonts}',   # LaTeX preamble, txfonts: nature font
        'grid.linestyle': '--', 'grid.linewidth': 0, 
        "xtick.direction":"in", 
        "xtick.major.width":0.8*alw, "xtick.major.size":1.0*tck_s, "xtick.minor.width":0.5*alw, "xtick.minor.size":tck_s*0.5, 
        "ytick.direction":"in", 
        "ytick.major.width":1.2*alw, "ytick.major.size":1.5*tck_s, "ytick.minor.width":0.6*alw, "ytick.minor.size":tck_s*0.8, 
        "axes.linewidth":alw
        }
    )

def plot_y3(): # y3 plot
    fig = plt.figure(figsize = (3.6*figs, 3.0*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lx, labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(ly, labelpad = lpad[1])

    xm = [-0.0, 10.0, 2.5]
    ax1.set_xlim(xm[0], xm[1]) # Range of axis
    ax1.set_xticks( np.arange(0.0, xm[1] + 1.e-3, xm[2]) ) # Ticks in xm[2] increments from xm[0] to xm[1]
    ax1.tick_params(axis = 'both', pad = tpad)

    ax1.set_yscale('log') ### log scale y axis

    ym = [np.nanmin(y3)*0.5, np.nanmax(y3)*2.0, None]
    ax1.set_ylim(ym[0], ym[1])
    int_logym = [math.ceil( np.log10(ym[0])*0.999 ), int( np.log10(ym[1]) )]
    y_tl = []; y_t = np.empty(0)
    for i in range(int_logym[0], int_logym[1] + 1): # Set y tick label
        y_t = np.append(y_t, 10**i)
        if i == 0:
            y_tl.append( r'$1$' )  # Set 10^0 as 1
        elif i == 1:
            y_tl.append( r'$10$' ) # Set 10^1 as 10
        else:
            y_tl.append( r'$10^{%d}$'%(i) )
    ax1.set_yticks(y_t, y_tl)

    ax1.axhline(1.0, lw = lw1*0.8, ls = 'dotted', label = r'$y = 1.0$') # Horizontal line
    xc = 4.8
    ax1.axvline(xc, lw = lw1*0.8, ls = 'dashed', dashes = [2, 5], color = color2, label = r'$x = x_\mathrm{c}$') # Vertical line, dashes: setting of dash line
    ax1.text((xc - xm[0])/(xm[1] - xm[0]) + 0.02 , 0.08, r'$x = x_\mathrm{c}$', color = color2, transform = ax1.transAxes) # Put text in graph
    ytext = y3[0]
    ax1.axhline(ytext, lw = lw1*0.8, c = 'C4', ls = 'dashdot')
    ax1.text(1.03, (np.log(ytext) - np.log(ym[0]))/(np.log(ym[1]) - np.log(ym[0])), r'$y_0 = %4.1f$'%(ytext), transform = ax1.transAxes)

    ### Plot loaded data
    ax1.plot(x1, y3, lw = lw1, ls = 'solid', marker = 'o', ms = ms1*1.2, mew = lw1*2., mfc = color0, color = color0, 
            alpha = 0.8,     # Transparency
            clip_on = False, # Plot out of graph with clip_on = False
            zorder = 8,      # Plot in front with larger zorder
            label = 'Normal font y3') 

    ### Legend setting
    # h1, l1 = ax1.get_legend_handles_labels()
    # ax1.legend(h1, l1, 
    # bbox_to_anchor = (1.05, 1.0), loc = "upper left", # bbox_to_anchor: position of loc
    # framealpha = 1.0, fancybox = False, 
    # edgecolor = "black").get_frame().set_linewidth(alw*0.8)

    ### Save as file
    fig.savefig(plotdir + "03plot_log" + ext, bbox_inches = "tight") # Remove margin with bbox_inches = "tight"

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)
    plot_y3()
    if fshow == 1:
        plt.show()
    print("end main")
