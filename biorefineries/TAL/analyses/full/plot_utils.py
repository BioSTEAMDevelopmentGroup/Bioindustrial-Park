#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

#%%
from biosteam.plots import plot_kde
from matplotlib import pyplot as plt
from matplotlib.colors import hex2color, LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator, LinearLocator, FixedLocator

#%%
def plot_kde_formatted(
                        xdata,
                        ydata,
                        
                        xlabel = r"$\bfMPSP$" + " " + r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
                        
                        ylabel = r"$\mathrm{\bfGWP}_{\bf100}$" + " " + r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$",
                           
                        xticks=None,
                        yticks=None,
                        
                        n_minor_ticks = 3,
                        
                        fig_width = 5,
                        fig_height = 5,
                        
                        show_x_ticks = True,
                        xbox_kwargs=dict(light=hex2color("#A97802"), dark=(0,0,0)),
                        ybox_kwargs=dict(light=hex2color('#607429'), dark=(0,0,0)),
                        
                        save_fig = True,
                        filename = 'Bivariate_KDE.png',
                        
                        vlines=[],
                        hlines=[],
                        
                        fill_between=[], # [(x_iterable, y1, y2, color), ...]
                        fill_between_zorders=[1000],
                        fill_between_alphas=[0],
                        plot_zorders=None,
                        cmaps=('Greens',),
                        
                        ax = None,
                        
                        ):
    axis_given = not (ax is None)
    fig, ax = plot_kde(xdata, 
                ydata,
                ax=ax,
                xticks=xticks,
                xticklabels=xticks,
                yticks=yticks,
                cmaps=cmaps,
                yticklabels=yticks,
                xbox_kwargs=xbox_kwargs,
                ybox_kwargs=ybox_kwargs,
                zorders=plot_zorders,
                )
    
    if not axis_given:
        if vlines:
            ax.vlines(vlines, yticks[0], yticks[-1], 
                      colors='#c0c1c2', 
                      linestyles='solid', zorder=0)
        if hlines:
            ax.vlines(hlines, xticks[0], xticks[-1],
                      colors='#c0c1c2', 
                      linestyles='solid', zorder=0)
        if fill_between:
            for fb, fbzo, fba in zip(fill_between, fill_between_zorders, fill_between_alphas):
                ax.fill_between(fb[0], fb[1], fb[2], 
                                # color='#c0c1c2',
                                # color='gray',
                                color=fb[3],
                                zorder=fbzo,
                                alpha=fba,
                                )
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        
        
        ####################### METHOD 1 ############################
        ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
        ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
        # ########--########
        ax.tick_params(
            axis='y',          # changes apply to the y-axis
            which='both',      # both major and minor ticks are affected
            direction='inout',
            right=True,
            width=0.65,
            )
        
        # ax.tick_params(
        #     axis='y',          # changes apply to the y-axis
        #     which='both',      # both major and minor ticks are affected
        #     direction='inout',
        #     right=False,
        #     width=0.65,
        #     )
        
        ax.tick_params(
            axis='y',          
            which='major',      
            length=7,
            )
        
        ax.tick_params(
            axis='y',          
            which='minor',      
            length=3.5,
            )
        
        ax.tick_params(
            axis='y',          
            which='major', 
            right=True,     
            length=7,
            )
        
        ax.tick_params(
            axis='y',          
            which='minor',     
            right=True,
            length=3.5,
            )
        
        
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            direction='inout',
            right=True,
            top=True,
            width=0.65,
            )
        
        # ax.tick_params(
        #     axis='x',          # changes apply to the x-axis
        #     which='both',      # both major and minor ticks are affected
        #     direction='inout',
        #     right=False,
        #     top=False,
        #     width=0.65,
        #     )
        
        ax.tick_params(
            axis='x',          
            which='major',      
            length=7,
            )
        
        ax.tick_params(
            axis='x',          
            which='minor',      
            length=3.5,
            )
        
        ax.tick_params(
            axis='x',          
            which='major',      
            length=7,
            right=True,
            top=True,
            )
        
        ax.tick_params(
            axis='x',          
            which='minor',      
            length=3.5,
            right=True,
            top=True,
            )
    ##################### END METHOD 1 ###################
    
    
    # #################### METHOD 2 #########################
    # ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
    
    
    
    # ax.tick_params(
    #     axis='y',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     direction='inout',
    #     # right=True,
    #     width=1,
    #     )
    
    # ax.tick_params(
    #     axis='y',          
    #     which='major',      
    #     length=5,
    #     )
    
    # ax.tick_params(
    #     axis='y',          
    #     which='minor',      
    #     length=3,
    #     )
    
    # ax2 = ax.twinx()
    
    # ax.set_ylim(min(yticks), max(yticks))
    # ax2.set_ylim(min(yticks), max(yticks))
    
    # if not yticks==[]:
    #     ax.set_yticks(yticks)
    #     ax2.set_yticks(yticks)
    #     l = ax.get_ylim()
    #     l2 = ax2.get_ylim()
    #     f = lambda x : l2[0]+(x-l[0])/(l[1]-l[0])*(l2[1]-l2[0])
    #     ticks = f(ax.get_yticks())
    #     ax2.yaxis.set_major_locator(FixedLocator(ticks))
    #     ax2.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
    
    # else:
    #     ax2.set_yticks(ax.get_yticks())
    #     l = ax.get_ylim()
    #     l2 = ax2.get_ylim()
    #     f = lambda x : l2[0]+(x-l[0])/(l[1]-l[0])*(l2[1]-l2[0])
    #     ticks = f(ax.get_yticks())
    #     ax2.yaxis.set_major_locator(FixedLocator(ticks))
    #     ax2.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
        
    #     # loc = LinearLocator(numticks = len(y_ticks))
    #     # ax2.set_yticks(y_ticks)
    #     # ax.yaxis.set_major_locator(loc)
    #     # ax2.yaxis.set_major_locator(loc)
    #     # nticks = len(y_ticks)
    #     # ax.yaxis.set_major_locator(LinearLocator(nticks))
    #     # ax2.yaxis.set_major_locator(loc)
        
    #     # ax2.set_yticks(np.linspace(ax2.get_yticks()[0], ax2.get_yticks()[-1], len(ax.get_yticks())))
        
    # ax2.tick_params(
    #     axis='y',          
    #     which='both',      
    #     direction='in',
    #     # right=True,
    #     labelright=False,
    #     width=1,
    #     )
    
    # ax2.tick_params(
    #     axis='y',          
    #     which='major',      
    #     length=5,
    #     right=True,
    #     )
    
    # ax2.tick_params(
    #     axis='y',          
    #     which='minor',      
    #     length=3,
    #     right=True,
    #     )
    
    # ################### END METHOD 2 ######################
    # ax.tick_params(axis='both', which='minor', zorder=2000)
    ax.set_axisbelow(False)
    ax.figure.set_size_inches(w=fig_width, h=fig_height)
    
    if save_fig:
        plt.savefig(filename, dpi=600, bbox_inches='tight',
                    facecolor=ax.get_facecolor(),
                    transparent=False,
                    )
    return ax

#%%

def get_cmap_from_colors(colors, N_levels=90, name='new_cmap'):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    """
    return LinearSegmentedColormap.from_list(name, colors, N_levels)

#%%
GG_colors_extended = {
    'green': '#7BBD84', 
    'orange':       '#E58835', 
    'yellow':      '#F7C652', 
    'blue':        '#63C6CE', 
    'pink':        '#F8858A', 
    'gray':        '#94948C', 
    'purple':        '#734A8C', 
    'lilac':        '#D1C0E1', 
    'darker blue':        '#648496', 
    'darker green': '#607429',
    'magenta': '#A100A1',
    'lighter magenta': '#FEC1FE'
    }

# GG_colors_extended = {k:hex2color(v) for k,v in GG_colors_extended.items()}