# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 08:02:27 2023
@author: lavanya
"""
import biosteam as bst
import thermosteam as tmo
from biosteam import Unit
from biosteam.units.decorators import cost
from biosteam.units.design_tools.cost_index import CEPCI_by_year
from biorefineries.oleochemicals.chemicals_baseline import chems

#Ref: Rules of Thumb
#Molar mass cut off = 200 = 0.2 kDa. The usual range is 0.01–1 kDa. 
#Pressure: (between UF and RO) = 0.3–1.4 MPa.
@cost(
      ID = 'Nanofilter',
      basis = 'Total filtration area',
      units='m^2', 
      cost=240*2*2.75, #here, 250$/m^2 is the cost of the membrane,
      #2 is the factor for skid mounted nanofiltration unit. 
      #2.75 factor includes labor and material costs excluding those related to instrumentation
      CE=1000, #CEPCI of 1000, as provided in book
      lb = 18, #lowest area possible
      ub = 22,# highest area possible
      n=1.0, #exponential factor
      S=1.0, #base size at which cost is provided (above)
      lifetime = 5 #typical lifetime of desalination NF membranes
      )

         
class Nanofiltration_oilcane_unit(bst.Unit,isabstract = True):
      _N_ins = 1
      _N_outs = 2
      _units = {
                'Total filtration area': 'm^2',
                }
      def __init__(self, ID='', ins=(), outs=(), thermo=None,
                   rejection_factors = None, #dict {'chemical_name',rejection_fraction}
                   volume_reduction_factor = None #provide a number to account for reduction in amount of water
                   ):
          Unit.__init__(self, ID, ins, outs, thermo)
          self.rejection_factors = rejection_factors
          self.volume_reduction_factor = volume_reduction_factor
          self.feed_rate_L_per_sec = 0.5 #based on L/s.m^2 provided in the book
         
            
      def _run(self):
        feed = self.ins[0]
        retenate,permeate, = self.outs
        if len(self.rejection_factors) == len(feed.available_chemicals) - 1:
            for i in self.rejection_factors:
                    mass_water = feed.imass['Water']
                    temp_chem_mass = feed.imass[i]
                    retenate.imass[i] = temp_chem_mass*(self.rejection_factors[i])
                    permeate.imass[i] = temp_chem_mass*(1-self.rejection_factors[i])
                    retenate.imass['Water'] = r_water_mass=  mass_water/self.volume_reduction_factor
                    permeate.imass['Water'] = mass_water - r_water_mass
        else:
            print('Please provide all the rejection factors except water')
                    
      def _design(self):
            m3_per_h_to_L_per_s = 0.2777
            process_vol_rate = self.ins[0].F_vol*m3_per_h_to_L_per_s#provides value in L/s
            self.design_results['Total filtration area'] = process_vol_rate/ self.feed_rate_L_per_sec
            
            
            
            
            
            
            
            
            
            
            
            