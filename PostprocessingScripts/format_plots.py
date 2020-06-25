import xarray as xr
import numpy as np
from numpy.linalg import inv, norm, eigh
from scipy.sparse import diags, identity
from scipy.stats import linregress
import pandas as pd
from tqdm import tqdm
import copy
from os.path import join

# clustering
from sklearn.cluster import KMeans

import math
import itertools

# Plotting
import matplotlib.pyplot as plt
from matplotlib import rcParams, colorbar, colors
from matplotlib.colors import LinearSegmentedColormap
import mpl_toolkits.mplot3d
from matplotlib.collections import PolyCollection, LineCollection
import cartopy.crs as ccrs
import cartopy
import cartopy.feature
from cartopy.mpl.patch import geos_to_path


# Figure sizes
SCALE = 2
# BASE_FIG_SIZE = 6
BASE_WIDTH = 8
BASE_HEIGHT = 4.5

# Fontsizes
TITLE_FONTSIZE = 18
SUBTITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10

# Position
TITLE_LOC = 1.1
CBAR_PAD = 0.05
LABEL_PAD = 20
CBAR_LABEL_PAD = 75

# Other font details
rcParams['font.family'] = 'serif'
rcParams['font.size'] = LABEL_FONTSIZE*SCALE
rcParams['text.usetex'] = True

def color(k, cmap='inferno', lut=10):
    c = plt.cm.get_cmap(cmap, lut=lut)
    return colors.to_hex(c(k))

def cmap_trans(cmap, ncolors=300, nalpha=20):
    color_array = plt.get_cmap(cmap)(range(ncolors))

    # change alpha values
    color_array[:,-1] = np.append(np.linspace(0.0, 1.0, nalpha),
                                  np.ones(ncolors-nalpha))

    # create a colormap object
    map_object = LinearSegmentedColormap.from_list(name=str(cmap) + '_trans',colors=color_array)

    return map_object

def cmap_trans_center(cmap, ncolors=300, nalpha=20):
    color_array = plt.get_cmap(cmap)(range(ncolors))

    # change alpha values
    half_l = math.floor((ncolors - nalpha)/2)
    half_r = math.ceil((ncolors - nalpha)/2)
    color_array[:,-1] = np.concatenate((np.ones(half_l),
                                        np.linspace(1, 0, nalpha/2),
                                        np.linspace(0, 1, nalpha/2),
                                        np.ones(half_r)))

    # create a colormap object
    map_object = LinearSegmentedColormap.from_list(name=str(cmap) + '_ctrans' ,
                                                   colors=color_array)

    return map_object

def get_figsize(aspect, rows, cols,
                max_width=BASE_WIDTH*SCALE,
                max_height=BASE_HEIGHT*SCALE):
    # if rows > 1 or cols > 1:
    #     factor = max(rows, cols)/2
    # else:
    #     factor = 1

    # if cols > 1:
    #     added_width = 1.25
    # else:
    #     added_width = 1

    # aspect *= added_width

    if aspect > 1: # width > height
        figsize = (max_width,
                   max_width/aspect)
    else: # width < height
        figsize = (max_height*aspect,
                   max_height)
    return figsize

def get_aspect(rows, cols, aspect=None, maps=False, lats=None, lons=None):
    if maps:
        aspect = np.cos(np.mean([np.min(lats), np.max(lats)])*np.pi/180)
        xsize = np.ptp([np.max(lons), np.min(lons)])*aspect
        ysize = np.ptp([np.max(lats), np.min(lats)])
        aspect = xsize/ysize
    return aspect*cols/rows

def make_axes(rows=1, cols=1, aspect=None, maps=False, lats=None, lons=None):
    aspect = get_aspect(rows, cols, aspect, maps, lats, lons)
    figsize = get_figsize(aspect, rows, cols)
    kw = {}
    if maps:
        kw['subplot_kw'] = {'projection' : ccrs.PlateCarree()}
    # if (rows + cols) > 2:
    #     kw['constrained_layout'] = True
        # figsize = tuple(f*1.5 for f in figsize)
    fig, ax = plt.subplots(rows, cols, figsize=figsize, **kw)
    # plt.subplots_adjust(right=1)
    return fig, ax

def add_cax(fig, ax):
    try:
        axis = ax[-1, -1]
        height = ax[0, -1].get_position().y1 - ax[-1, -1].get_position().y0
    except IndexError:
        axis = ax[-1]
        # height = ax[-1].get_position().height
        height = ax[0].get_position().y1 - ax[-1].get_position().y0
    except TypeError:
        axis = ax
        height = ax.get_position().height

    # x0
    cbar_pad_inches = 0.75
    fig_width = fig.get_size_inches()[0]
    x0_init = axis.get_position().x1
    x0 = (fig_width*x0_init + cbar_pad_inches)/fig_width

    # y0
    y0 = axis.get_position().y0

    # Width
    cbar_width_inches = 0.2
    width = cbar_width_inches/fig_width

    # Make axis
    cax = fig.add_axes([x0, y0, width, height])

    return cax

def get_figax(rows=1, cols=1, aspect=1,
              maps=False, lats=None, lons=None, kw={}):
    n = len(kw)
    if 'figax' in kw.keys():
        fig, ax = kw.pop('figax')
    else:
        fig, ax = make_axes(rows, cols, aspect, maps, lats, lons)

    if (rows > 1) or (cols > 1):
        for axis in ax.flatten():
            axis.set_facecolor('0.98')
        # plt.subplots_adjust(hspace=0.1, wspace=0.4)
    else:
        ax.set_facecolor('0.98')

    if n > 0:
        return fig, ax, kw
    else:
        return fig, ax

def add_labels(ax, xlabel, ylabel,
               fontsize=LABEL_FONTSIZE*SCALE,
               labelsize=TICK_FONTSIZE*SCALE,
               labelpad=LABEL_PAD):
    ax.set_xlabel(xlabel, fontsize=fontsize, labelpad=labelpad)
    ax.set_ylabel(ylabel, fontsize=fontsize, labelpad=labelpad)
    ax.tick_params(axis='both', which='both', labelsize=labelsize)
    return ax

def add_legend(ax, frameon=False, fontsize=LABEL_FONTSIZE*SCALE, **kw):
    ax.legend(frameon=frameon, fontsize=fontsize, **kw)
    return ax

def add_title(ax, title, y=TITLE_LOC, fontsize=TITLE_FONTSIZE*SCALE, **kw):
    ax.set_title(title, y=y, fontsize=fontsize, **kw)
    return ax

def get_square_limits(xdata, ydata):
    # Get data limits
    dmin = min(np.min(xdata), np.min(ydata))
    dmax = max(np.max(xdata), np.max(ydata))
    pad = (dmax - dmin)*0.05
    dmin -= pad
    dmax += pad

    try:
        # get lims
        ylim = kw.pop('ylim')
        xlim = kw.pop('xlim')
        xy = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
    except:
        # set lims
        xlim = ylim = xy = (dmin, dmax)

    return xlim, ylim, xy, dmin, dmax

def format_map(ax, lats, lons, fontsize=TICK_FONTSIZE*SCALE, draw_labels=True):
    ax.set_ylim(min(lats), max(lats))
    ax.set_xlim(min(lons), max(lons))
    ax.add_feature(cartopy.feature.OCEAN, facecolor='0.98')
    ax.add_feature(cartopy.feature.LAND, facecolor='0.98')
    ax.coastlines(color='grey')
    gl = ax.gridlines(linestyle=':', draw_labels=draw_labels, color='grey')
    gl.xlabel_style = {'fontsize' : fontsize}
    gl.ylabel_style = {'fontsize' : fontsize}
    return ax

def format_cbar(cbar, cbar_title=''):
    # cbar.set_label(cbar_title, fontsize=BASEFONT*SCALE,
    #                labelpad=CBAR_LABEL_PAD)
    cbar.ax.text(5.25, 0.5, cbar_title,
                 ha='center', va='center',
                 rotation='vertical',
                 fontsize=LABEL_FONTSIZE*SCALE,
                 transform=cbar.ax.transAxes)
    cbar.ax.tick_params(axis='both', which='both',
                        labelsize=TICK_FONTSIZE*SCALE)
    return cbar

def plot_one_to_one(ax):
    xlim, ylim, _, _, _ = get_square_limits(ax.get_xlim(),
                                            ax.get_ylim())
    ax.plot(xlim, xlim, c='0.1', lw=2, ls=':',
            alpha=0.5, zorder=0)
    return ax

def save_fig(fig, loc, name):
    fig.savefig(join(loc, name + '.png'), bbox_inches='tight', dpi=300)
    print('Saved %s' % name + '.png')
