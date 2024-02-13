"""
Animation, LaTeX font
12/02/2024
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import seaborn as sns
import os
import sys
import warnings
warnings.simplefilter("ignore", category = DeprecationWarning)
import time
start = time.time()

fshow = 0   # Display popup 0: no, 1: yes
### Output setting
fsize = 1   # 0: small size for paper, 1: large size for talk
ext = 'pdf' # File extension (pdf, svg, png, etc.)
ext_ani = 'mov'  # Animation extension (mov, mp4, etc.)
plotdir = 'plot' # Plot file save directry

datadir = 'data'
ddir_ani = '04animation'
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
ddir_ani, ext_ani = datadir + ddir_ani + '/', '.' + ext_ani

### Data for plot
nt_p = 51 # number of time step
tmax = 10.0
dt_p = tmax/(nt_p - 1) # time per step
nxi  = np.array([128, 64])    # Number of x and y data
Lxi  = np.array([4.0, 2.0]) # Max of x and y
dxi  = Lxi/nxi              # Grid length
### Create data. Can comment out if data files exist
# os.makedirs(ddir_ani, exist_ok = True) # Make data animation dir
# dfile = ddir_ani + 'data_xy.dat'
# f2 = open(dfile, 'w')
# for it in range(nt_p):
#     dfile   = ddir_ani + 'data%5.5d.dat'%(it)
#     f = open(dfile, 'w')
#     for i2 in range(nxi[1] + 1):
#         y = dxi[1]*i2
#         for i1 in range(nxi[0] + 1):
#             x = dxi[0]*i1
#             cxy = np.sin(2*math.pi*(y/Lxi[1] - x/Lxi[1] + it*dt_p/3.0))
#             f.write('%8.3f\n'%(cxy))
#             if it == 0:
#                 f2.write('%8.3f %8.3f\n'%(x, y))
#     f.close()
# f2.close()

if os.path.exists(datadir) == False:  # Cancel if datadir or plotdir does not exist
    sys.exit("datadir error: %s"%datadir)
os.makedirs(plotdir, exist_ok = True) # Make plotdir
if os.path.exists(plotdir) == False:
    sys.exit("plotdir error: %s"%plotdir)

### Load text file
x1  = np.loadtxt(ddir_ani + 'data_xy.dat',       usecols = 0, dtype = 'float32') # usecols: column number, dtype (float32: single precision, float64: double precision, int64: integer, etc.)
y1  = np.loadtxt(ddir_ani + 'data_xy.dat',       usecols = 1, dtype = 'float32')
cxy0= np.loadtxt(ddir_ani + 'data%5.5d.dat'%(0), usecols = 0, dtype = 'float32')
x1  =  x1.reshape([nxi[1] + 1, -1]) # Reshape loaded data into 2D array. y and x grid points. -1: automatically determined
y1  =  y1.reshape([nxi[1] + 1, -1]) # Reshape loaded data into 2D array. y and x grid points. -1: automatically determined

### Setting of size, label, etc.
lx, ly, lc = r'$x/L_0$', r'$y/L_0$', r'$\tilde{c}_\mathrm{a}(x,y,t)$' # Use TeX character with r''
lc2 = r'$\alpha(x,y,t)$'
if fsize == 0:
    figs, figs_ani = 1.0, 0.3
    fs1 = 1.0
    tck_s1, alw = 3, 0.625
    ctxt1 = 'paper'
    lpad = 5 # Space between tick label and label
    tpad = 5 # Space between axis and tick label
else:
    figs, figs_ani = 2.0, 1.0
    fs1 = 1.5
    tck_s1, alw = 7, 1.25
    ctxt1, ext = 'talk', '_talk' + ext
    lpad = 11
    tpad = 11

vm1 = [-1.0, 0.6] # Range of color bar
vm1.append( (vm1[1] - vm1[0])/4.0 )
vm2 = [math.sqrt(10.0)*1.0e-2, math.sqrt(10.0)*1.0e2] # Range of color bar

alpha_max = 1.0e2
def alpha_t(array1, t_in):
    return (array1 - np.nanmin(cxy0))*alpha_max*math.cos(math.pi*t_in/5.0)**2.0 + 8.0e-2

hide_ani_cmap = 1 ### Must be 1
if hide_ani_cmap == 1: ### Setting for color map animation
    it_m, cxy_m = [0, 0], [np.nanmin(cxy0), np.nanmax(cxy0)]
    extend_ani2 = 'neither'
    for it in range(nt_p):
        cxy  = np.loadtxt(ddir_ani + 'data%5.5d.dat'%(it), usecols = 0, dtype = 'float32')
        if np.nanmin(cxy) < cxy_m[0]:
            it_m[0], cxy_m[0] = it, np.nanmin(cxy)
        if np.nanmax(cxy) > cxy_m[1]:
            it_m[1], cxy_m[1] = it, np.nanmax(cxy)

        cxy2 = np.loadtxt(ddir_ani + 'data%5.5d.dat'%(it), usecols = 0, dtype = 'float32')
        cxy2 = alpha_t(cxy2, dt_p*it)
        if extend_ani2 == 'neither':
            if vm2[0] > np.nanmin(cxy2):
                extend_ani2 = 'min'
            if vm2[1] < np.nanmax(cxy2):
                extend_ani2 = 'max'
        if ( extend_ani2 == 'min' and vm2[1] < np.nanmax(cxy2) ) or ( extend_ani2 == 'max' and vm2[0] > np.nanmin(cxy2) ):
            extend_ani2 = 'both'

    extend_ani = 'neither' # extend setting. If data exceed min or max value, make edge triangle
    if vm1[0] > cxy_m[0] and vm1[1] < cxy_m[1]:
        extend_ani = 'both'
    elif vm1[0] > cxy_m[0]:
        extend_ani = 'min'
    elif vm1[1] < cxy_m[1]:
        extend_ani = 'max'

    int_logvm = [math.ceil( np.log10(vm2[0])*0.999 ), int( np.log10(vm2[1]) )]
    c2_tl = []; c2_t = np.empty(0)
    for i in range(int_logvm[0], int_logvm[1] + 1): # Set y tick label
        c2_t = np.append(c2_t, 10**i)
        if i == 0:
            c2_tl.append( r'$1$' )  # Set 10^0 as 1
        elif i == 1:
            c2_tl.append( r'$10$' ) # Set 10^1 as 10
        else:
            c2_tl.append( r'$10^{%d}$'%(i) )

    xm_c, ym_c = [0, Lxi[0], 1.0], [0, Lxi[1], 1.0]
    extent1 = [xm_c[0] - dxi[0]*0.5, xm_c[1] + dxi[0]*0.5, ym_c[0] - dxi[1]*0.5, ym_c[1] + dxi[1]*0.5] # Range of color map. dxi: grid length of xi. Careful for imshow format

hide_define_cmap = 1 ### Must be 1
if hide_define_cmap == 1: ### Define and create cmap
    ### Existing cmap
    cmap1 = 'viridis'
    ### Create original cmap
    #### Intermittent change
    cnmax = 10
    cim = cm.get_cmap('Blues', cnmax) # Create data 'Blues' divided into cnmax
    clist = [cim(0)]
    for i in range(1, cim.N):
        clist.append(cim(i))
    mycmap1 = ListedColormap(clist)
    #### Irregular change
    cnmax = 401
    cim = cm.get_cmap('viridis', cnmax) # Create data 'viridis' divided into cnmax
    cnum, tp = [0, 100, 200, 300, cim.N - 1], [50, 150, 250, 350] # cnum: number of colors to use, tp: Where to switch colors
    sm = 10 # Smoothness of color change (transition length around tp[j])
    clist, j = [ cim(cnum[0]) ], 0
    for i in range(1, cim.N):
        cn0 = 0
        if j < len(tp):
            if i >= tp[j] + sm:
                j = j + 1
            elif abs( i - tp[j] ) <= sm:
                cn0 = int( (cnum[j + 1] - cnum[j])/(2*sm)*(i - (tp[j] - sm)) )
        clist.append( cim(cnum[j] + cn0) )
    mycmap2 = LinearSegmentedColormap.from_list('colormap_name', clist)

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
        "xtick.direction":"in", "xtick.major.width":0.8*alw, "xtick.major.size":tck_s, 
        "ytick.direction":"in", "ytick.major.width":0.8*alw, "ytick.major.size":tck_s, 
        "axes.linewidth":alw
        }
    )

class mk_animation:
    def update_single_cmap(self, i, ax1, fig):
        global cbar_initialized
        t_p = dt_p*i + 1.0e-6
        ax1.cla()
        ax1.set_aspect('equal') # Follow aspect ratio
        ax1.set_xlabel(lx, labelpad = lpad) # Axis label
        ax1.set_ylabel(ly, labelpad = lpad)
        ax1.tick_params(axis = 'both', pad = tpad)
        ax1.set_xlim(xm_c[0], xm_c[1]) # Range of axis
        ax1.set_ylim(ym_c[0], ym_c[1])
        ax1.set_xticks( np.arange(xm_c[0], xm_c[1] + 1.0e-3, xm_c[2]) ) # x ticks in xm[2] increments from xm[0] to xm[1]
        ax1.set_yticks( np.arange(ym_c[0], ym_c[1] + 1.0e-3, ym_c[2]) ) # y ticks
        ax1.text(0.45, 1.07, r'$t = %4.1f \,\mathrm{s}$'%(t_p), verticalalignment = 'center_baseline', transform = ax1.transAxes)

        c_plt = np.loadtxt(ddir_ani + 'data%5.5d.dat'%(i), usecols = 0, dtype = 'float32')
        c_plt = c_plt.reshape([nxi[1] + 1, -1]) # Reshape loaded data into 2D array. y and x grid points. -1: automatically determined
        im = ax1.imshow(c_plt,                     # c_plt is 2D array
                        interpolation = 'bicubic', # Interpolation (bilinear, none, etc.)
                        extent = extent1, 
                        cmap = mycmap1,            # Color map type
                        origin = 'lower',          # Set origin lower
                        norm = Normalize(vmin = vm1[0], vmax = vm1[1])
                        )

        if not cbar_initialized:
            cbar_initialized = True ### Important for color map animation! 
            axpos = ax1.get_position() # Get ax1 position. x0: left, x1: right, y0: bottom, y1: top, height: height
            pp_ax = fig.add_axes([axpos.x1 + 0.02, axpos.y0, 0.03, axpos.height]) # Color bar left bottom x and y, width, height
            pp = fig.colorbar(im, ax = ax1, orientation = "vertical", cax = pp_ax, extend = extend_ani)
            pp.ax.yaxis.set_tick_params(pad = tpad, right = False, which = "minor")
            pp.set_ticks( np.arange(vm1[0], vm1[1] + 1.e-3, vm1[2]) ) # Color bar ticks
            pp.set_label(lc, labelpad = lpad, loc = 'center', rotation = 90)  # Color bar label, rotate with rotation (Default: 90)

    def ani_single_cmap(self): ### Single color map
        fig = plt.figure(figsize = (16.0*figs_ani, 9.0*figs_ani), dpi = 100, linewidth = 0)
        ax1 = fig.add_subplot(111)
        chartB = ax1.get_position()
        ax1.set_position([chartB.x0, chartB.y0, chartB.width*0.88, chartB.height]) # Graph position and size

        ani_t = tmax*1.5 # Animation time
        ani = animation.FuncAnimation(fig, self.update_single_cmap, fargs = (ax1, fig), interval = ani_t*1000/(nt_p - 1), frames = nt_p)
        ani.save(plotdir + '04plot_ani_single' + ext_ani, writer = 'ffmpeg')

    def update_multiple(self, i, ax0, ax_l, axc, fig):
        global cbar_initialized
        t_p = dt_p*i + 1.0e-8
        ax0.cla()
        ax0.axis("off")
        ax0.text(0.1, 1.0, r'$t = %4.1f \,\mathrm{s}$'%(t_p), verticalalignment = 'center_baseline', transform = ax0.transAxes)
        for i1 in range(len(axc)):
            axc[i1].cla()
            axc[i1].set_aspect('equal') # Follow aspect ratio
            axc[i1].set_xlabel(lx, labelpad = lpad) # Axis label
            axc[i1].set_ylabel(ly, labelpad = lpad)
            axc[i1].tick_params(axis = 'both', pad = tpad)
            axc[i1].set_xlim(xm_c[0], xm_c[1]) # Range of axis
            axc[i1].set_ylim(ym_c[0], ym_c[1])
            axc[i1].set_xticks( np.arange(xm_c[0], xm_c[1] + 1.0e-3, xm_c[2]) ) # x ticks in xm[2] increments from xm[0] to xm[1]
            axc[i1].set_yticks( np.arange(ym_c[0], ym_c[1] + 1.0e-3, ym_c[2]) ) # y ticks

            if i1 == 0:
                c_plt = np.loadtxt(ddir_ani + 'data%5.5d.dat'%(i), usecols = 0, dtype = 'float32')
                norm1 = Normalize(vmin = vm1[0], vmax = vm1[1])
                cm_ml = mycmap1
                p_arr1 = np.arange(vm1[0], vm1[1] + 1.e-5, vm1[2])
                pp_list = [extend_ani, p_arr1, np.round(p_arr1, 3), lc]
            else:
                c_plt = np.loadtxt(ddir_ani + 'data%5.5d.dat'%(i), usecols = 0, dtype = 'float32')
                c_plt = alpha_t(c_plt, t_p)
                norm1 = LogNorm(vmin = vm2[0], vmax = vm2[1]) # Log scale
                cm_ml = mycmap2
                pp_list = [extend_ani2, c2_t, c2_tl, lc2]

            c_plt = c_plt.reshape([nxi[1] + 1, -1])    # Reshape loaded data into 2D array. y and x grid points. -1: automatically determined
            im = axc[i1].imshow(c_plt,                     # c_plt is 2D array
                            interpolation = 'bicubic', # Interpolation (bilinear, none, etc.)
                            extent = extent1, 
                            cmap = cm_ml,            # Color map type
                            origin = 'lower',          # Set origin lower
                            norm = norm1
                            )

            if not cbar_initialized:
                if i1 == len(axc) - 1:
                    cbar_initialized = True ### Important for color map animation! 
                axpos = axc[i1].get_position() # Get ax1 position. x0: left, x1: right, y0: bottom, y1: top, height: height
                pp_ax = fig.add_axes([axpos.x1 + 0.015, axpos.y0, 0.015, axpos.height]) # Color bar left bottom x and y, width, height
                pp = fig.colorbar(im, ax = axc[i1], orientation = "vertical", cax = pp_ax, extend = pp_list[0])
                pp.ax.yaxis.set_tick_params(pad = tpad, right = False, which = "minor")
                pp.ax.tick_params(axis = 'y', length = tck_s1*0.8, width = alw, direction = 'inout')
                pp.set_ticks(pp_list[1]) # Color bar ticks
                pp.set_ticklabels(pp_list[2]) # Color bar tick labels
                pp.set_label(pp_list[3], labelpad = lpad, loc = 'center', rotation = 90)  # Color bar label, rotate with rotation (Default: 90)

            ax_l[i1].cla()
            if i1 == 0:
                ax_l[i1].set_xlabel(lx, labelpad = lpad) # Axis label
                ax_l[i1].set_ylabel(r'$\tilde{c}_\mathrm{a} (x,y=0,t)$', labelpad = lpad)
                ax_l[i1].set_xlim(xm_c[0], xm_c[1])
                ax_l[i1].set_ylim(-1.1,    1.1)
                ax_l[i1].set_xticks( np.arange(xm_c[0], xm_c[1] + 1.0e-3, xm_c[2]) )
                ax_l[i1].plot(x1[0, :], c_plt[0, :], lw = 1.5*fs1, c = 'black')
            else:
                ax_l[i1].set_xlabel(ly, labelpad = lpad)
                ax_l[i1].set_ylabel(r'$\alpha(x=0,y,t)$', labelpad = lpad)
                ax_l[i1].set_xlim(ym_c[0], ym_c[1])
                ax_l[i1].set_xticks( np.arange(ym_c[0], ym_c[1] + 1.0e-3, ym_c[2]) )
                ax_l[i1].set_yscale('log') ### log scale y axis
                ax_l[i1].set_ylim(vm2[0],  vm2[1])
                ax_l[i1].set_yticks(c2_t, c2_tl)
                ax_l[i1].plot(y1[:, 0], c_plt[:, 0], lw = 1.5*fs1)
            ax_l[i1].tick_params(axis = 'both', pad = tpad)

    def ani_multiple(self): ### Multiple
        fig = plt.figure(figsize = (16.0*figs_ani, 9.0*figs_ani), dpi = 100, linewidth = 0)
        ax0 = fig.add_subplot(111)

        axl_1 = fig.add_subplot(2, 2, 2)
        axl_2 = fig.add_subplot(2, 2, 4)
        ax_l = [axl_1, axl_2]
        ax_l_s = 0.75
        for i in range(len(ax_l)):
            chartB = ax_l[i].get_position()
            ax_l[i].set_position([chartB.x0 + 0.10, chartB.y0 + chartB.height*(1.0 - ax_l_s)/2.0, chartB.width*ax_l_s, chartB.height*ax_l_s]) # Graph position and size

        axc_1 = fig.add_subplot(2, 2, 1)
        axc_2 = fig.add_subplot(2, 2, 3)
        axc = [axc_1, axc_2]
        for i in range(len(axc)):
            chartB = axc[i].get_position()
            axc[i].set_position([chartB.x0 - 0.03, chartB.y0, chartB.width*0.85, chartB.height]) # Graph position and size

        ani_t = tmax*2.5 # Animation time
        ani = animation.FuncAnimation(fig, self.update_multiple, fargs = (ax0, ax_l, axc, fig), interval = ani_t*1000/(nt_p - 1), frames = nt_p)
        ani.save(plotdir + '04plot_ani_multiple' + ext_ani, writer = 'ffmpeg')

##=================== main ===================##
if __name__ == '__main__':
    minute = int((time.time() - start)/60)
    print(minute, 'm%6.2fs, start main '%(time.time() - start - minute*60.0))
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    print('c_plt min', it_m[0], cxy_m[0], ', max', it_m[1], cxy_m[1])
    sns_set(fs1, tck_s1, alw, ctxt1)

    cbar_initialized = False
    mk_animation().ani_single_cmap()
    minute = int((time.time() - start)/60)
    print(minute, 'm%6.2fs, end animation 1 '%(time.time() - start - minute*60.0))

    cbar_initialized = False
    mk_animation().ani_multiple()

    if fshow == 1:
        plt.show()
    minute = int((time.time() - start)/60)
    print(minute, 'm%6.2fs, end main '%(time.time() - start - minute*60.0))
