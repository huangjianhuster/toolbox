import numpy as np
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


def contour_plot(dat, bins=200, colormap=None, color="#00BFFF", contour_levels=None, save_path=None):
    """
    to plot 2d histogram distribution with contours
        by default, 2d histogram density will be calculated and use the "color" parameter to highlight the most dense reion.
        colors will decay linearly from the "color" parameter to "white", i.e., from "#00BFFF" decays to "white"
    dat: shape of (2, dim_y)
    color: the color for the most dense region; the #00BFFF is the color skyblue
    return None
    """
    if not colormap:
        # By default, use linear decay color scheme
        custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', [(0, 'white'), (1, color)])
    else:
        custom_cmap = colormap

    # processing raw data
    heatmap, xedges, yedges = np.histogram2d(dat[0,:], dat[1,:], bins=[bins, bins], density=True)
    extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]

    # plot
    fig, axs = plt.subplots(1, 1, figsize=(6, 6)) # ,  subplot_kw=dict(aspect='equal'))
    im = axs.imshow(heatmap.T, origin='lower', extent=extent, cmap=custom_cmap, norm=mpl.colors.LogNorm(),\
                    zorder=2, alpha=0.8)

    # set xlim and ylim
    # axs.set_xlim(x_min, x_max)
    # axs.set_ylim(y_min, y_max)

    cax = fig.add_axes([axs.get_position().x1+0.01, axs.get_position().y0, 0.02, axs.get_position().height])
    colorbar = fig.colorbar(im, cax=cax)
    # colorbar.set_label('Density') #, labelpad=20)
    colorbar.ax.tick_params(labelsize=8)

    # add contour
    # mask the heatmap since there are lots of zeros
    heatmap = np.ma.masked_where(heatmap == 0, heatmap)
    if not contour_levels:
        contour_levels = np.logspace(np.log10(np.min(heatmap)), np.log10(np.max(heatmap)), 5)
    else:
        contour_levels = contour_levels
    print("contour_levels: ", contour_levels)

    # exponential decay of the linewidth
    linewidths = np.logspace(np.log10(0.1), np.log10(3), len(contour_levels))
    contours = axs.contour(heatmap.T, extent=extent, levels=contour_levels, \
                            colors=color, linewidths=linewidths, zorder=5)

    # Add contour levels as lines on colorbar
    # cb_ticks = colorbar.ax.get_yticks()
    # for tick,thickiness in zip(contour_levels, linewidths):
    #     colorbar.ax.axhline(tick, color=color, linestyle='-', linewidth=thickiness)
    #     colorbar.ax.text(3.5, tick, f'{tick:.1e}', color='b', ha='right', va='center', fontsize=8)
    # axs.grid(zorder=0)
    # plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()


def ramachandran_contour(dat, bins=200, color="#00BFFF", contour_levels=None, save_path=None):
    """
    to plot ramachandran distribution
        by default, 2d histogram density will be calculated and use the "color" parameter to highlight the most dense reion.
        colors will decay linearly from the "color" parameter to "white", i.e., from "#00BFFF" decays to "white"
    dat: shape of (2, dim_y)
    color: the color for the most dense region; the #00BFFF is the color skyblue
    return None
    """
    # Linear decay color scheme
    custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', [(0, 'white'), (1, color)])
    
    # processing raw data
    heatmap, xedges, yedges = np.histogram2d(dat[0,:], dat[1,:], bins=[bins, bins], density=True)

    # plot
    fig, axs = plt.subplots(1, 1, figsize=(6, 5),  subplot_kw=dict(aspect='equal'))
    im = axs.imshow(heatmap.T, origin='lower', extent=[xedges[0], xedges[-1],\
                    yedges[0], yedges[-1]], cmap=custom_cmap, norm=mpl.colors.LogNorm(),\
                    zorder=2, alpha=0.8)
    colorbar = fig.colorbar(im, ax=axs, fraction=0.046, pad=0.02)
    colorbar.set_label('Density', labelpad=20)
    
    # add contour
    # mask the heatmap since there are lots of zeros
    heatmap = np.ma.masked_where(heatmap == 0, heatmap)
    if not contour_levels:
        contour_levels = np.logspace(np.log10(np.min(heatmap)), np.log10(np.max(heatmap)), 5)
    else:
        contour_levels = contour_levels
    print("contour_levels: ", contour_levels)

    # exponential decay of the linewidth
    linewidths = np.logspace(np.log10(0.1), np.log10(3), len(contour_levels))
    contours = axs.contour(heatmap.T, extent=[xedges[0], xedges[-1],\
                        yedges[0], yedges[-1]], levels=contour_levels, \
                            colors=color, linewidths=linewidths, zorder=5)

    # Add contour levels as lines on colorbar
    cb_ticks = colorbar.ax.get_yticks()
    for tick,thickiness in zip(contour_levels, linewidths):
        colorbar.ax.axhline(tick, color=color, linestyle='-', linewidth=thickiness)
        # colorbar.ax.text(3.5, tick, f'{tick:.1e}', color='b', ha='right', va='center', fontsize=8)

    axs.set_xlabel(r'$\phi$')
    axs.set_ylabel(r'$\psi$')
    axs.set_xlim([-180, 180])
    axs.set_ylim([-180, 180])
    axs.set_xticks(np.arange(-180, 240, 60))
    axs.set_yticks(np.arange(-180, 240, 60))
    axs.grid(zorder=0)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()
