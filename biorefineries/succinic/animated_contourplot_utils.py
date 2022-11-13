# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 10:15:45 2022

@author: sarangbhagwat
"""

#%% Plot MPSP

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import imageio


def generate_animated_contourplot(w_data_vs_x_y_at_multiple_z, # shape = z * x * y
                                  x_data,
                                  y_data,
                                  z_data,
                                  x_label, # title of the x axis
                                  y_label, # title of the y axis
                                  z_label, # title of the z axis
                                  w_label, # title of the color axis
                                  x_ticks,
                                  y_ticks,
                                  z_ticks,
                                  w_ticks, # labeled, lined contours (a subset of w_levels)
                                  w_levels, # unlabeled, filled contour areas (labeled and ticked only on color bar)
                                  x_units,
                                  y_units,
                                  z_units,
                                  w_units,
                                  w_tick_width=0.5, # width for labeled, lined contours
                                  fmt_clabel = lambda cvalue: "{:.2f}".format(cvalue), # format of contour labels
                                  gridspec_kw={'height_ratios': [1, 20]},
                                  fontname={'fontname':'Arial'},
                                  figwidth=4.5,
                                  dpi=600,
                                  cmap='viridis',
                                  z_marker_color='b',
                                  z_marker_type='v',
                                  axis_title_fonts={'size': {'x': 12, 'y':12, 'z':12, 'w':12},},
                                  gap_between_figures=20., 
                                  fps=3, # animation frames (z values traversed) per second
                                  n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                  animated_contourplot_filename='animated_contourplot'
                                 ):
    results = np.array(w_data_vs_x_y_at_multiple_z)
    
    def create_frame(z_index):
        fig, axs = plt.subplots(2, 1, constrained_layout=True, 
                                gridspec_kw=gridspec_kw)
        fig.set_figwidth(figwidth)
        ax = axs[0]
        a = [z_data[z_index]]
        ax.hlines(1,1,1)
        ax.set_xlim(min(z_ticks), max(z_ticks))
        ax.set_ylim(0.5,1.5)
        ax.xaxis.tick_top()
        ax.set_xlabel(z_label + " [" + z_units + "]",  
                      fontsize=axis_title_fonts['size']['z'],
                      **fontname)
        ax.set_xticks(z_ticks,
                      **fontname)
    
        y = np.ones(np.shape(a))
        ax.plot(a,y,
                color=z_marker_color, 
                marker=z_marker_type,
                ms = 7,)
        ax.axes.get_yaxis().set_visible(False)
        
        ax = axs[1]
        im = ax.contourf(x_data, y_data, results[z_index],
                          cmap=cmap,
                         levels=w_levels,
                         )
  
        clines = ax.contour(x_data, y_data, results[z_index],
                   levels=w_ticks,
                    colors='black',
                   linewidths=w_tick_width)
        
        ax.clabel(clines, 
                   w_ticks,
                   fmt=fmt_clabel, 
                  fontsize=10,
                  colors='black',
                  )
        
        ax.set_ylabel(y_label + " [" + y_units + "]",  
                      fontsize=axis_title_fonts['size']['y'],
                      **fontname)
        ax.set_xlabel(x_label + " [" + x_units + "]", 
                      fontsize=axis_title_fonts['size']['x'],
                      **fontname)
        
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        
        cbar = plt.colorbar(im, cax=cax, 
                     ticks = w_levels).set_label(label=w_label + " [" + w_units + "]", 
                                                           size=axis_title_fonts['size']['w'],
                                                           loc='center',
                                                           **fontname
                                                           )
        ax.set_title(' ', fontsize=gap_between_figures)

        plt.savefig(f'./TRY_images/MPSP_TRY_prod_{z_data[z_index]}.png', 
                    transparent = False,  
                    facecolor = 'white',
                    bbox_inches='tight',
                    dpi=dpi,
                    )                                
        plt.close()
        
        
    for z_index in range(len(z_data)):
        create_frame(z_index)
          
    frames = []
    for z_index in range(len(z_data)):
        image = imageio.v2.imread(f'./TRY_images/MPSP_TRY_prod_{z_data[z_index]}.png')
        frames.append(image)
    
    
    if n_loops==('inf' or 'infinite' or 'infinity' or np.inf):
        imageio.mimsave('./' + animated_contourplot_filename + '.gif',
                        frames,
                        fps=fps,
                        ) 
    else:
        imageio.mimsave('./' + animated_contourplot_filename + '.gif',
                        frames,
                        fps=fps,
                        loop=n_loops,
                        ) 
        
    
