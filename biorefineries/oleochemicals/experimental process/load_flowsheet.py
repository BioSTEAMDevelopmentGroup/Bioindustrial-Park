"""
Created on Tue Jul 26 15:32:26 2022

@author: LENOVO
"""

from biorefineries import oleochemicals as oc
from importlib import reload
reload(oc)
oc.load('azelaic_acid')
oc.sys.diagram(number = True)
oc.sys.simulate()