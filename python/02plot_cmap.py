"""
Color map etc. LaTeX font
10/01/2023
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import os
import sys

fshow = 0 # 結果のポップアップ表示
### 入出力の設定
fsize = 0   # 0：小さいサイズ（論文用）　1：大きいサイズ（プレゼン用）
ext = 'pdf' # 保存ファイルの拡張子　pdf,svg,pngなど
plotdir = 'plot' # 保存用ディレクトリ
os.makedirs(plotdir, exist_ok = True) # ディレクトリ作成

datadir = 'data'
plotdir, datadir, ext = plotdir + '/', datadir + '/', '.' + ext
dfile   = datadir + 'data2.dat'

if os.path.exists(plotdir) == False:  # plotdirやdatadirが存在しないときに中止する
    sys.exit("plotdir error: %s"%plotdir)
elif os.path.exists(datadir) == False:
    sys.exit("datadir error: %s"%datadir)

### プロット用のデータ
nxi = np.array([16, 8]) # データの数
Lxi = np.array([4.0, 2.0]) # x, yの最大値
dxi = Lxi/nxi
### プロット用データの作成。すでにファイルが存在しているときはコメントアウトしてOK
# f = open(dfile, 'w')
# for i2 in range(nxi[1] + 1):
#     y = dxi[1]*i2
#     for i1 in range(nxi[0] + 1):
#         x = dxi[0]*i1
#         cxy = np.sin(2*math.pi*(y/Lxi[1] - x/Lxi[1]))
#         f.write('%8.3f %8.3f %8.3f\n'%(x, y, cxy))
# f.close()

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
x1  = np.loadtxt(dfile, usecols = 0, dtype = 'float64') # usecolsは列番号　dtypeは実数float64, 整数int64など
y1  = np.loadtxt(dfile, usecols = 1, dtype = 'float64')
cxy = np.loadtxt(dfile, usecols = 2, dtype = 'float64')
#### 2次元形式に変換
x1  =  x1.reshape([nxi[1] + 1, -1]) # 読み込んだデータを2次元配列に変換。y, xの要素数。-1で自動
y1  =  y1.reshape([nxi[1] + 1, -1]) # 読み込んだデータを2次元配列に変換。y, xの要素数。-1で自動
cxy = cxy.reshape([nxi[1] + 1, -1]) # 読み込んだデータを2次元配列に変換。y, xの要素数。-1で自動

### サイズ、ラベルなどの設定
lx, ly, lc = r'$x$', r'$y$', r'$c_\mathrm{a}(x,y)$' # r''でTeX文字にできる
if fsize == 0:
    figs = 1.
    fs1 = 1.
    tck_s1, alw = 3, 0.625
    ctxt1 = 'paper'
    lpad = [5, 5] # 軸とラベルの間隔
    tpad = [3, 5] # 軸と数値の間隔
else:
    figs = 2.
    fs1 = 1.5
    tck_s1, alw = 7, 1.25
    ctxt1, ext = 'talk', '_talk' + ext
    lpad = [10, 8]
    tpad = [8, 12]

### 既存のcmap
cmap1 = 'viridis'
### cmapの自作
#### 断続的な変化
cnmax = 10
cim = cm.get_cmap('Blues', cnmax) # 'Blues'をcnmax分割したデータを取得
clist = [cim(0)]
for i in range(1, cim.N):
    clist.append(cim(i))
mycmap1 = ListedColormap(clist)
#### 不規則に変化させる
cnmax = 401
cim = cm.get_cmap('viridis', cnmax) # 'viridis'をcnmax分割したデータを取得
cnum, tp = [0, 100, 200, 300, cim.N - 1], [50, 150, 250, 350] # cnum使う色の番号、tp色を切り替える場所の指定
sm = 10 # 変化の滑らかさ（0からcnumの間隔の半分くらいまでの値）
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


def plot_cxy(): # cxyプロット用
    fig = plt.figure(figsize = (4.8*figs, 2.4*figs), dpi = 100, linewidth = 0)
    ax1 = fig.add_subplot(111)
    ax1.set_aspect('equal') # アスペクト比を数値通りに表示する
    ax1.spines["top"].set_linewidth(alw)
    ax1.spines["left"].set_linewidth(alw)
    ax1.spines["bottom"].set_linewidth(alw)
    ax1.spines["right"].set_linewidth(alw)

    ax1.set_xlabel(lx, labelpad = lpad[0]) # 軸ラベル
    ax1.set_ylabel(ly, labelpad = lpad[1])
    xm, ym = [0, Lxi[0], 1], [0, Lxi[1], 1]
    ax1.set_xlim(xm[0], xm[1]) # 軸の範囲
    ax1.set_ylim(ym[0], ym[1])
    ax1.set_xticks( np.arange(xm[0], xm[1] + 1.e-3, xm[2]) ) # xm[0]からxm[1]までxm[2]刻みの目盛り線
    ax1.set_yticks( np.arange(ym[0], ym[1] + 1.e-3, ym[2]) ) # yの目盛り線
    # ax1.set_xticklabels([0, 1, 2, 3, r'$x_\mathrm{max}$']) # 目盛りに文字列を入れる
    ax1.tick_params(axis='x', pad = tpad[0])
    ax1.tick_params(axis='y', pad = tpad[1])

    c_plt = cxy
    # c_plt = x1 # デバック用（x1, y1でプロットしてimshowの軸の方向など確認）
    # c_plt = y1
    # print(c_plt)
    vm = [np.nanmin(c_plt), np.nanmax(c_plt), (np.nanmax(c_plt) - np.nanmin(c_plt))/4.] # カラーバーの範囲
    print('vmin, vmax:', vm)
    extent1 = [xm[0] - dxi[0]*0.5, xm[1] + dxi[0]*0.5, ym[0] - dxi[1]*0.5, ym[1] + dxi[1]*0.5] # グラフの範囲の指定。dxi xiの刻み幅。imshowの表示方法に注意する
    ### 読み込みデータのプロット
    im = ax1.imshow(c_plt, # cxyは2次元配列
                    interpolation = 'bicubic', # 補間方法 bilinear, noneなど
                    extent = extent1, 
                    cmap = cmap1, # カラーマップの種類 'viridis'など
                    origin = 'lower', # 原点を下にする
                    vmin = vm[0], vmax = vm[1])

    extend1 = 'neither' # extendの設定。データの値がカラーバーの値の範囲を超えていれば三角にする
    if vm[0] > np.nanmin(c_plt) and vm[1] < np.nanmax(c_plt):
        extend1 = 'both'
    elif vm[0] > np.nanmin(c_plt):
        extend1 = 'min'
    elif vm[1] < np.nanmax(c_plt):
        extend1 = 'max'
    axpos = ax1.get_position() # グラフの位置情報を取得　x0左端, x1右端, y0下端, y1上端, height高さ
    pp_ax = fig.add_axes([axpos.x1 + 0.02, axpos.y0, 0.03, axpos.height]) # カラーバーの左下のx,y,幅,高さ
    pp = fig.colorbar(im, ax = ax1, orientation = "vertical", cax = pp_ax, extend = extend1)
    cticks = np.arange(vm[0], vm[1] + 1.e-3, vm[2])
    pp.set_ticks(cticks) # カラーバーの目盛り線
    # pp.set_ticklabels([cticks[0], cticks[1], cticks[2], r'$c_\mathrm{c}$', cticks[4]]) # 目盛りに文字列を入れる
    pp.set_label(lc, labelpad = 5, loc = 'center', rotation = 90)  # カラーバーのラベル。rotationでラベルを回転（デフォルトは90）

    fig.savefig(plotdir + "02plot_cxy" + ext, bbox_inches = "tight") # 保存。bbox_inches="tight"で余白をなくす

##=================== main ===================##
if __name__ == '__main__':
    print("start main")
    print("datadir: ", datadir, ", plotdir: ", plotdir, ", ctxt: ", ctxt1)
    sns_set(fs1, tck_s1, alw, ctxt1)
    plot_cxy()
    if fshow == 1:
        plt.show()
    print("end main")
