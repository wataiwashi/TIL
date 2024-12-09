"""
Lorenz sync
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

### Lorenz parameter
sigma = 10.0
b = 8.0/3.0
r = 28.0
# r = 0.8
dt = 1.0e-5
tmax, dt_plt = 50.0, 2.0e-5
tmin = 30.0

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
        dot_xyz[1] = r*xyz[0] - xyz[1] - 20.0*xyz[0]*xyz[2]
        dot_xyz[2] = 5.0*xyz[0]*xyz[1] - b*xyz[2]
        xyz[:] = xyz[:] + dot_xyz[:]*dt
    return t_plt, xyz_plt

def lorenz_receive(x0, y0, z0, u_drive, tmin, tmax, dt_plt):
    n_p = int((tmax + 1.0e-8)/dt_plt) + 1 - int(tmin/dt_plt + 1.0e-6)
    t_plt = np.zeros(n_p)
    xyz_plt = np.zeros((3, n_p))
    xyz = np.array([x0, y0, z0])
    dot_xyz = np.zeros(3)
    for nt in range(n_p):
        t = nt*dt_plt
        i_p = int((t + 1.0e-8)/dt_plt)
        t_plt[i_p] = t + tmin
        xyz_plt[:, i_p] = xyz[:]

        dot_xyz[0] = -sigma*xyz[0] + sigma*xyz[1]
        dot_xyz[1] = r*u_drive[nt] - xyz[1] - 20.0*xyz[0]*xyz[2] # input u_drive
        dot_xyz[2] = 5.0*xyz[0]*xyz[1] - b*xyz[2]
        xyz[:] = xyz[:] + dot_xyz[:]*dt_plt
    return t_plt, xyz_plt

t_plt1, xyz_plt1 = lorenz(0.0, 1.0, 0.0, tmax, dt_plt)
t_plt1, xyz_plt1 = t_plt1[int(tmin/dt_plt + 1.0e-6):], xyz_plt1[:, int(tmin/dt_plt + 1.0e-6):]
t_plt2, xyz_plt2 = lorenz_receive(0.0, 0.01, 0.0, xyz_plt1[0], tmin, tmax, dt_plt)

hide_0 = 1 ### Must be 1
if hide_0 == 1: ### Setting of size, label, etc. 
    lt, lx, ly, lz = r'$t$', r'$X$', r'$v, v_r$', r'$Z$' # Use TeX character with r''
    if fsize == 0:
        figs = 1.
        fs1, lw1, ms1 = 1., 1.2, 6.
        tck_s1, alw = 3, 0.625
        ctxt1 = 'paper'
        lpad = [5, 5] # Space between tick label and label
        tpad = [5, 5] # Space between axis and tick label
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

def plot_y():
    fig = plt.figure(figsize = (4.8*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lt, labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(ly, labelpad = lpad[1])
    xm = [tmin, tmax]
    ax1.set_xlim(xm[0], xm[1]) # Range of axis
    ax1.set_ylim(-4, 4)
    ax1.tick_params(axis = 'x', pad = tpad[0])
    ax1.tick_params(axis = 'y', pad = tpad[1])

    ax1.plot(t_plt1, xyz_plt1[1], lw = lw1, ls = 'solid', color = 'C0', alpha = 1.0, clip_on = False, zorder = 8, label = r'$v(t=%4.1f$'%tmin + r'$)=%5.2f$'%xyz_plt1[1, 0])

    ax1.plot(t_plt2, xyz_plt2[1], lw = lw1*1.5, ls = 'solid', color = 'C1', alpha = 0.7, clip_on = False, zorder = 8, label = r'$v_r(t=%4.1f$'%tmin + r'$)=%5.2f$'%xyz_plt2[1, 0])

    ax1.scatter(t_plt1[0], xyz_plt1[1, 0], marker = 'o', s = 40, c = 'C0', clip_on = False, zorder = 10)
    ax1.scatter(t_plt2[0], xyz_plt2[1, 0], marker = 'o', s = 40, c = 'C1', clip_on = False, zorder = 10)

    ax1.set_xticks(ax1.get_xticks(), set_ticks_digit(ax1.get_xticks()))
    ax1.set_yticks(ax1.get_yticks(), set_ticks_digit(ax1.get_yticks()))

    ### Legend
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ## Save file
    fig.savefig(plotdir + "lorenz_sync" + ext, bbox_inches = "tight")

def plot_error():
    fig = plt.figure(figsize = (4.8*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lt, labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(r'$|\bm{e}|$', labelpad = lpad[1])
    xm = [tmin, tmax]
    ax1.set_xlim(xm[0], xm[1]) # Range of axis
    ax1.set_yscale('log')
    ax1.tick_params(axis = 'x', pad = tpad[0])
    ax1.tick_params(axis = 'y', pad = tpad[1])

    error = np.sqrt((xyz_plt1[0, :] - xyz_plt2[0, :])**2 + (xyz_plt1[1, :] - xyz_plt2[1, :])**2 + (xyz_plt1[2, :] - xyz_plt2[2, :])**2)
    ax1.plot(t_plt1, error, lw = lw1, ls = 'solid', color = 'C0', alpha = 1.0, clip_on = False, zorder = 8, label = r'$v(t=30)=%6.3f$'%xyz_plt1[1, 0])

    ax1.set_xticks(ax1.get_xticks(), set_ticks_digit(ax1.get_xticks()))
    arr_yt, lis_y = set_ticks_digit_log(ax1.get_ylim())
    ax1.set_yticks(arr_yt, lis_y)

    ## Save file
    fig.savefig(plotdir + "lorenz_sync_error" + ext, bbox_inches = "tight")

def plot_message():
    fig = plt.figure(figsize = (4.8*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lt, labelpad = lpad[0]) # Axis label
    ax1.set_ylabel(r'$m, m_r$', labelpad = lpad[1])
    xm = [tmin, tmax]
    ax1.set_xlim(xm[0], xm[1]) # Range of axis
    ax1.tick_params(axis = 'x', pad = tpad[0])
    ax1.tick_params(axis = 'y', pad = tpad[1])

    message = np.zeros(len(t_plt1))
    t_sin, f, A = t_plt1[0], 0.0, 0.0
    for i in range(1, len(message)):
        if t_plt1[i] >= t_sin:
            T = random.uniform(0.5, 1.5)
            t_sin = T + t_plt1[i]
            f = random.randint(1, 4)
            A = random.uniform(0.0, 0.002)
        message[i] = A*np.sin((t_plt1[i] - t_sin)*float(f)*2.0*np.pi/T)

    u_drive = xyz_plt1[0] + message
    t_plt, xyz_plt = lorenz_receive(u_drive[0], xyz_plt1[1, 0], xyz_plt1[2, 0], u_drive, tmin, tmax, dt_plt)

    ax1.plot(t_plt1, message, lw = lw1, ls = 'solid', color = 'C0', alpha = 1.0, zorder = 8, label = r'$m$')

    ax1.set_ylim(ax1.get_yticks()[0], ax1.get_yticks()[-1])

    ax1.plot(t_plt1, u_drive - xyz_plt[0], lw = lw1*1.5, ls = 'solid', color = 'C1', alpha = 0.7, zorder = 8, label = r'$m_r$')

    ax1.set_xticks(ax1.get_xticks(), set_ticks_digit(ax1.get_xticks()))
    ax1.set_yticks(ax1.get_yticks(), set_ticks_digit(ax1.get_yticks()))

    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    fig.savefig(plotdir + "message" + ext, bbox_inches = "tight")

def plot_xyz():
    fig = plt.figure(figsize = (3.6*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111, projection = '3d')
    ax1.set_xlabel(lx) # Axis label
    ax1.set_ylabel(ly)
    ax1.set_zlabel(lz)

    ax1.scatter(xyz_plt1[0, 0], xyz_plt1[1, 0], xyz_plt1[2, 0], marker = '*', c = 'red')
    ax1.plot(xyz_plt1[0], xyz_plt1[1], xyz_plt1[2], lw = lw1*0.25, ls = 'solid', color = 'C0', alpha = 1.0, clip_on = False, zorder = 8, label = r'$y_0=%10.5f$'%xyz_plt1[1, 0])
    # ax1.plot(xyz_plt2[0], xyz_plt2[1], xyz_plt2[2], lw = lw1*0.25, ls = 'solid', color = 'C1', alpha = 0.7, clip_on = False, zorder = 8, label = r'$y_0=%10.5f$'%xyz_plt2[1, 0])

    ### Legend
    # h1, l1 = ax1.get_legend_handles_labels()
    # ax1.legend(h1, l1, bbox_to_anchor = (1.0, 1.0), loc = "upper left", framealpha = 1.0, fancybox=False, edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ### save
    fig.savefig(plotdir + 'lorenz.jpg', dpi = 600, bbox_inches = 'tight', pad_inches = 0.3)

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("plotdir: ", plotdir, ", ctxt: ", ctxt1)
    print('ntmax, n_p: ', int((tmax + 1.0e-8)/dt), int((tmax + 1.0e-8)/dt_plt))
    print('tmin: ', tmin, ', x0, y0, z0: ', np.round(xyz_plt1[:, 0], 4))
    sns_set(fs1, tck_s1, alw, ctxt1)

    plot_y()
    plot_error()
    plot_message()
    # plot_xyz()

    if fshow == 1:
        plt.show()
    print("end main")
