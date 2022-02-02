# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 15:44:35 2022

@author: saran
"""

from biorefineries.TAL.system_TAL import *
from matplotlib import pyplot as plt
from pandas import DataFrame


def get_MPSP_for_feedstock_flow(flow_rate):
    feedstock.F_mass = flow_rate
    try:
        return get_SA_MPSP(), get_feedstock_dry_mass(), SA.F_mass, R302.ins[1].imass['Glucose', 'Xylose'].sum()
    except:
        print('Error at feedstock flow rate of %s kg/h; returning 0.'%str(flow_rate))
        return 0., 0., 0.
    
def get_feedstock_dry_mass():
    return feedstock.F_mass - feedstock.imass['Water']

#%% Run
feedstock_flow_rates = np.linspace(10000, sum(baseline_feedflow)*5., 100)
metrics_data = [get_MPSP_for_feedstock_flow(i) for i in feedstock_flow_rates]

#%%% Manage data
metrics_data_dict = {'TAL MPSP': [i[0] for i in metrics_data],
                     'Feedstock dry mass': [i[1]*24./1000. for i in metrics_data],
                     'TAL capacity': [i[2]*24./1000. for i in metrics_data],
                     'Glucose and xylose to fermentation': [i[3]*24./1000. for i in metrics_data]
                     }

metrics_data_df = DataFrame(metrics_data_dict, feedstock_flow_rates*24./1000.)

metrics_data_df.to_excel('TAL_MPSP_vs_capacity.xlsx')
#%% Plot MPSP vs Feedstock dry mass
x = 'Feedstock dry mass'
y = 'TAL MPSP'
fig, ax = plt.subplots()
ax.plot(metrics_data_dict[x], metrics_data_dict[y])
# ax = metrics_data_df.plot(xlabel='Incidence rate [%]', ylabel='Test predictive value [%]',
#                  )
ax.set_xlabel(x + ' [metric tonne/d]')
ax.set_ylabel(y + ' [$/kg]')
ax.set_ylim(ymin = 0)
ax.set_xlim(xmin = 0)
plt.show()

#%% Plot MPSP vs TAL capacity
x = 'TAL capacity'
y = 'TAL MPSP'
fig, ax = plt.subplots()
ax.plot(metrics_data_dict[x], metrics_data_dict[y])
# ax = metrics_data_df.plot(xlabel='Incidence rate [%]', ylabel='Test predictive value [%]',
#                  )
ax.set_xlabel(x + ' [metric tonne/d]')
ax.set_ylabel(y + ' [$/kg]')
ax.set_ylim(ymin = 0)
ax.set_xlim(xmin = 0)
plt.show()

#%% Plot MPSP vs Glucose and xylose to fermentation
x = 'Glucose and xylose to fermentation'
y = 'TAL MPSP'
fig, ax = plt.subplots()
ax.plot(metrics_data_dict[x], metrics_data_dict[y])
# ax = metrics_data_df.plot(xlabel='Incidence rate [%]', ylabel='Test predictive value [%]',
#                  )
ax.set_xlabel(x + ' [metric tonne/d]')
ax.set_ylabel(y + ' [$/kg]')
ax.set_ylim(ymin = 0)
ax.set_xlim(xmin = 0)
plt.show()


