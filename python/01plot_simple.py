"""
Line plot etc. LaTeX font
10/01/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

### 入出力の設定
fsize = 0   # 0：小さいサイズ（論文用）　1：大きいサイズ（プレゼン用）
ext = 'pdf' # 保存ファイルの拡張子　pdf,svg,pngなど
plotdir = 'plot' # 保存用ディレクトリ
os.makedirs(plotdir, exist_ok = True) # ディレクトリ作成

datadir = 'data'
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
dfile   = datadir + 'data1.dat'

if os.path.exists(plotdir) == False or os.path.exists(datadir) == False: # plotdirやdatadirが存在しないときに中止する
    sys.exit("Error")

def sns_set(fs, tck_s, alw, ctxt):
    sns.set(
        context = ctxt,  # フォントサイズ・線幅（'paper' or 'talk'）
        palette = sns.color_palette("colorblind"),
        font = "serif",       # フォント指定
        # font = "sans-serif",  # サンセリフ体
        font_scale = fs,  # フォントスケール指定（これを変えるとcontextで決まるプリセットを更にいじれる）
        style = None, 
        # style = "whitegrid",   # 背景白，グリッドあり
        rc = {"text.usetex": True, 
        'text.latex.preamble' : r'\newcommand{\ssmr}[1]{_\mathrm{#1}}', # LaTeXのプリアンブル
        # 'text.latex.preamble' : r'\usepackage{txfonts}',   # nature系のフォント、ギリシャ文字あるとき
        'grid.linestyle': '--', 'grid.linewidth': 0, 
        "xtick.direction":"in", "xtick.major.width":0.8*alw, "xtick.major.size":tck_s, 
        "ytick.direction":"in", "ytick.major.width":0.8*alw, "ytick.major.size":tck_s, 
        "axes.linewidth":alw
        }
    )

### 読み込み
x1 = np.loadtxt(dfile, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など
y1 = np.loadtxt(dfile, usecols = 1, dtype = 'float64')
y2 = np.loadtxt(dfile, usecols = 2, dtype = 'float64')
y3 = np.loadtxt(dfile, usecols = 3, dtype = 'float64')

### サイズ、ラベルなどの設定
lx, ly = r'Distance, $x$', r'$y$' # r''でTeX文字にできる
if fsize == 0:
    fs1, lw1, ms1 = 1., 1., 6.
    tck_s1, alw = 3, 0.625
    ctxt1 = 'paper'
    lpad = [5, 5] # 軸とラベルの間隔
    tpad = [3, 5] # 軸と数値の間隔
else:
    fs1, lw1, ms1 = 1.5, 3., 14.
    tck_s1, alw = 7, 1.25
    ctxt1, ext = 'talk', '_talk' + ext
    lpad = [10, 8]
    tpad = [8, 12]


def plot_y1y2(): # y1, y2プロット用
    if fsize == 0:
        fig = plt.figure(figsize = (3.6, 2.4), dpi = 100, linewidth = 0)
    else:
        fig = plt.figure(figsize = (7.2, 4.8), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    ax1.set_ylabel(ly, labelpad = lpad[1])
    xm, ym = [0, 10], [-35, 10]
    ax1.set_xlim(xm[0], xm[1]) # 軸の範囲
    ax1.set_ylim(ym[0], ym[1])
    ax1.set_xticks(np.arange(0, xm[1] + 1.e-3, 2.5)) # xmaxまで2.5刻みの目盛り線
    ax1.tick_params(axis='x', pad = tpad[0])
    ax1.tick_params(axis='y', pad = tpad[1])

    ### 水平線と鉛直線
    ax1.axhline(0, lw = lw1*0.8, ls = 'dotted', label = r'$y=0$')
    ax1.axvline(4.8, lw = lw1*0.8, ls = 'dashed', dashes = [2, 5], color = 'C2', label = r'$x=x\ssmr{c}$') # dashesで破線の間隔などを設定できる
    ytext = -22.
    ax1.axhline(ytext, lw = lw1*0.8, c = 'C4', ls = 'dashdot')
    ax1.text(0.7, (ytext - ym[0])/(ym[1] - ym[0]) + 0.05, r'$y=%4.1f$'%(ytext), transform = ax1.transAxes) # 図中にtextを入れる

    ### 読み込みデータのプロット
    ax1.plot(x1, y1, lw = lw1, ls = 'solid', marker = 'o', ms = ms1, color = 'C0', 
            alpha = 0.6,     # 透明度
            clip_on = False, # プロット枠外にもプロットする
            zorder = 8,      # zorderが大きいほど前面に表示される
            label = 'Normal font y1') 
    ax1.plot(x1, y2, lw = lw1, ls = 'solid', ms = ms1, marker = 's', color = 'C1', alpha = 1, clip_on = False, zorder = 7, label = r'$y_2$')

    ### 凡例の設定
    h1, l1 = ax1.get_legend_handles_labels()
    ax1.legend(h1, l1, 
    bbox_to_anchor = (1.0, 1.0), loc = "upper left", # bbox_to_anchorは凡例のlocの座標
    framealpha = 1.0, fancybox=False, 
    edgecolor = "black").get_frame().set_linewidth(alw*0.8)
    ### 保存
    fig.savefig(plotdir + "01plot_y1y2" + ext, bbox_inches = "tight") # bbox_inches="tight"で余白をなくす

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)
    plot_y1y2()
    print("end main")
