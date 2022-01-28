# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 08:17:38 2021
@author: yrc2
"""
import biosteam as bst
import thermosteam as tmo

class OzonolysisReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
    
    @property
    def effluent(self):
        return self.outs[0]
    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent
        
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=17, N=None, V=None, T=373.15, P=101325,
                 Nmin=2, Nmax=36):
        
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                   tau = tau , N = N, V = V, T = T, 
                                   P = P ,Nmin = Nmin , Nmax = Nmax)
        
#When you inherit, you write all the parameters of the base class. Why is Nmin = Nmin
#What is this - should it be ozochemicals

        c = self.chemicals      
        
        self.Oleic_acid_conversion = 0.808
        #Oleic_acid formula = C18H34O2
        #EPS formula = C18H34O3
        
        #self.nonanal_selectivity = 0.094 
        #mol EPS/ mol OA
                             
        self.selectivity_oxiraneoctanoic_acid = 0.071
        #self.Oxononanoic_acid_selectivity = 0.029
        #'C9H16O3'
        
        #self.Nonanoic_acid_selectivity = 0.277 
        #C9H18O2
        
        
        self.selectivity_Azelaic_acid = 0.442 * 2 
        #carbon moles of EPS*(1 mol of EPS/ 9 C moles)/ carbon moles of oleic acid*(1 mol of OA/ 18 C moles)
        # AA C9H16O4
        
    def _setup(self):
        # selectivity_Azelaic_acid = selectivity_Nonanoic_acid = X2 * X1
        # selectivity_Nonanal = selectivity_Oxonanoic_acid = X1 * (1 - X2)
        # selectivity_oxiraneoctanoic_acid,_3-octyl- = (1 - X1)
        
        #Considering that 
        X=self.Oleic_acid_conversion
        X1 = 1 - (self.selectivity_oxiraneoctanoic_acid)
        X2 = self.selectivity_Azelaic_acid / X1
        
        self.reactions = tmo.SeriesReaction([
            tmo.Rxn('Oleic_acid + H2O2 -> Epoxy_stearic_acid + Water ', 'Oleic_acid', X=self.Oleic_acid_conversion),
            tmo.Rxn('Epoxy_stearic_acid + H2O2 -> Nonanal + Oxononanoic_acid + H2O', 'oxiraneoctanoic_acid,_3-octyl-', X = X1),
            tmo.Rxn('Nonanal + Oxononanoic_acid + 2H2O2 -> Azelaic_acid + Nonanoic_acid+ 2H2O', 'Nonanal', X = X2)
            ])
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)
              
        self.reactions(effluent) 
        
        effluent.T = self.T
        effluent.P = self.P
        
# =============================================================================
# class D1(bst.BinaryDistillation):
#     
#     _N_ins = 1
#     _N_outs = 2   
#    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                 P= 3333.05921,   
#                 T= 503.15, LHK = ('Nonanoic_acid','Azelaic_acid'),
#                 k=2,
#                 Rmin=0.3,
#                 Lr=None,
#                 Hr=None,
#                 y_top = 0.8,
#                 x_bot= 0.2,                 
#                 product_specification_format='Composition',
#                 vessel_material='Carbon steel',
#                 tray_material='Carbon steel',
#                 tray_type='Sieve',
#                 tray_spacing=450,
#                 stage_efficiency=None,
#                 velocity_fraction=0.8,
#                 foaming_factor=1.0,
#                 open_tray_area_fraction=0.1,
#                 downcomer_area_fraction=None,
#                 is_divided=False,
#                 vacuum_system_preference='Liquid-ring pump',
#                 condenser_thermo=None,
#                 boiler_thermo=None,
#                 partial_condenser=True,
#                 
#                 ):
#         
#           self.P = P
#           self.T = T
#           self.product_specification_format = product_specification_format
#           self.y_top = y_top
#           self.x_bot = x_bot  
#                     
#           bst.BinaryDistillation.__init__(self, ID, ins, outs, thermo,
#                 P=P, LHK=LHK, k=k,
#                 Rmin=0.3,
#                 Lr=Lr,
#                 Hr=Hr,
#                 y_top= y_top,
#                 x_bot=x_bot, 
#                 product_specification_format=None,
#                 vessel_material='Carbon steel',
#                 tray_material='Carbon steel',
#                 tray_type='Sieve',
#                 tray_spacing=450,
#                 stage_efficiency=None,
#                 velocity_fraction=0.8,
#                 foaming_factor=1.0,
#                 open_tray_area_fraction=0.1,
#                 downcomer_area_fraction=None,
#                 is_divided=False,
#                 vacuum_system_preference='Liquid-ring pump',
#                 condenser_thermo=None,
#                 boiler_thermo=None,
#                 partial_condenser=True,
#         )
#                  
#     def _run(self):
#         feed = self.ins[0]
#         products = self.outs[0]
#         products.copy_like(feed)
#         products.T = self.T
#         products.P = self.P
# =============================================================================

# =============================================================================
# class Separator(bst.Unit):
#     #Why does this not have an _init_?
#     
#     _N_outs = 6
#         
#     def _run(self):
#         feed = self.ins[0]
#         IDs = ['Oleic_acid','Nonanal','Nonanoic_acid','Azelaic_acid', 
#                   'Oxononanoic_acid', 'oxiraneoctanoic_acid,_3-octyl-']
#         outs = self.outs[0]
#         
#         for stream, ID in zip(self.outs,IDs):
#             stream.imol[ID] = feed.imol[ID]
# =============================================================================

            

                   
                        
      
        
        
        
    
    
    
    
        