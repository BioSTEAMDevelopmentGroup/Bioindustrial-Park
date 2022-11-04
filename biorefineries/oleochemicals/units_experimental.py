"""
Created on Fri Oct 29 08:17:38 2021
@author: yrc2
"""
import biosteam as bst
import thermosteam as tmo
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream, MultiStream, equilibrium

class OxidativeCleavageReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 2
        
    def _setup(self):
        super()._setup() 
#TODO: change the reaction conversions
        Epoxide_formation = SRxn([Rxn('Oleic_acid + H2O2   -> Epoxy_stearic_acid ', 'Oleic_acid', X= 0.808),
                                  Rxn('Epoxy_stearic_acid + H2O -> DHSA', 'Epoxy_stearic_acid', X = 0.8),
                                  Rxn('DHSA + H2O2 -> Nonanal + Oxononanoic_acid ', 'DHSA', X = 0.8)])
        
        Side_reactions_1 =    SRxn([Rxn('Nonanal + H2O2 ->  Nonanoic_acid', 'Nonanal', X = 1),
                                    Rxn('Nonanoic_acid ->  Octane + CO2', 'Nonanoic_acid', X = 1),
                                    Rxn('Octane + H2O2 -> Octanal', 'Octane', X = 0.4)])
        Side_reactions_2 = PRxn([ Rxn('Oxononanoic_acid -> CO2 + Octanal', 'Oxononanoic_acid', X = 0.2),
                                  Rxn('Oxononanoic_acid + H2O2 ->  Azelaic_acid', 'Oxononanoic_acid', X = 0.8),
                                ])  
        Side_reactions_3 = PRxn([ Rxn('Octanal + H2O2 -> Octanoic_acid ', 'Octanal', X = 0.2),
                                  Rxn('Octanoic_acid + H2O2 -> Heptanal + CO2', 'Octanoic_acid',X = 0.4),
                                  Rxn('Heptanal + H2O2 -> Heptanoic_acid','Heptanal', X = 0.5)])
        
        # linoleic acid conversion to azelaic acid reported in: #https://doi.org/10.1016/j.indcrop.2022.115139
        Impurities_side_reactions = PRxn([Rxn('Linoleic_acid -> Azelaic_acid + Malonic_acid + Hexanoic_acid','Linoleic_acid', X = 0.8)
                                         ])  
    
        oxidative_cleavage_rxnsys = RxnSys(Epoxide_formation,
                                           Side_reactions_1,
                                           Side_reactions_2,
                                           Side_reactions_3,
                                           Impurities_side_reactions)
        self.reactions = oxidative_cleavage_rxnsys
        
    def _run(self):
        feed = self.ins[0]
        vent,effluent = self.outs
        effluent.copy_like(feed)
        self.reactions(effluent) 
        vent.copy_flow(effluent,'CO2',remove = True)
        effluent.copy_like(effluent)
        effluent.T = self.T
        effluent.P = self.P
        
class HeterogeneousReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 2
        
    def _setup(self):
        super()._setup() 
#Below reaction conversions are based on experimental results
#TODO: modify reactions below to change the reactant conversions
#Can add hydrogen peroxide decomposition reaction as well
        Epoxide_formation = SRxn([Rxn('Oleic_acid + H2O2   -> Epoxy_stearic_acid ', 'Oleic_acid', X= 0.684),
                                  Rxn('Epoxy_stearic_acid + H2O -> DHSA', 'Epoxy_stearic_acid', X = 0.8),
                                  Rxn('DHSA + H2O2 -> Nonanal + Oxononanoic_acid ', 'DHSA', X = 0.8)])
        
        Side_reactions_1 =    SRxn([Rxn('Nonanal + H2O2 ->  Nonanoic_acid', 'Nonanal', X = 1),
                                   ])
        Side_reactions_2 = PRxn([Rxn('Oxononanoic_acid + H2O2 ->  Azelaic_acid', 'Oxononanoic_acid', X = 0.8),
                                ])  

        oxidative_cleavage_rxnsys = RxnSys(Epoxide_formation,
                                           Side_reactions_1,
                                           Side_reactions_2,
                                           )
        self.reactions = oxidative_cleavage_rxnsys
        
    def _run(self):
        feed = self.ins[0]
        vent,effluent = self.outs
        effluent.copy_like(feed)
        self.reactions(effluent) 
        vent.copy_flow(effluent,'CO2',remove = True)
        effluent.copy_like(effluent)
        effluent.T = self.T
        effluent.P = self.P

    

class AACrystalliser(bst.units.BatchCrystallizer):
  
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,  
                 T = None
                  ):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                        tau=2, V=1e6, T=T)
        # https://www.alfa.com/en/catalog/B21568/
        # https://www.chembk.com/en/chem/Nonanoic%20acid#:~:text=Nonanoic%20acid%20-%20Nature%20Open%20Data%20Verified%20Data,yellow.%20Slightly%20special%20odor.%20Melting%20Point%2011-12.5%20%C2%B0c.
        self.AA_molefraction_330_15K = 0.0006996
        self.AA_molefraction_280_15K = 0.0000594
        self.NA_solubility_gL_at20DEGC = 0.267/1000
#Nonanoic acid melting point is 12.5
                        
    @property
    def Hnet(self):
        effluent = self.outs[0]
        solids = effluent['s']
        H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
        return H_out 

    def solubility(self, T):
        delta_T = 330.15 - 280.15
        delta_S = self.AA_molefraction_330_15K - self.AA_molefraction_280_15K
        m = delta_S/delta_T
        b = m*330.15
        c = self.AA_molefraction_330_15K - b
        S = m*T + c
        return S
    
#Assuming inlet at saturation. Therefore adding the feed at saturation of 330.15
    def _run(self):
        feed = self.ins[0]    
        outlet = self.outs[0]
        outlet.copy_like(feed)
        outlet.phases = ('s', 'l')
        x = self.solubility(self.T)
        outlet.sle('Azelaic_acid',
                    solubility=x,
                    T = self.T)
        outlet.imass['s','Nonanoic_acid'] = feed.imass['Nonanoic_acid']

        
    
# #Solubility data
# #Pelargonic acid is insoluble in water
# #Interpolation for Azelaic acid in water
# # m = 22-2/50-20 = 1.96
# # S = mT + C
# # S = 1.96T - 76.0
# #Therefore, we assume azelaic acid solubility at 95 deg cel is 110g/l

#TODO: Use the below while writing reaction conversions
        # self.reactions = tmo.SeriesReaction([
        #     tmo.Rxn('Oleic_acid + H2O2   -> Epoxy_stearic_acid + Water ', 'Oleic_acid', X=self.Oleic_acid_conversion),
        #     tmo.Rxn('Epoxy_stearic_acid + H2O2 -> Nonanal + Oxononanoic_acid + H2O', 'oxiraneoctanoic_acid,_3-octyl-', X = X1),
        #     tmo.Rxn('Nonanal + Oxononanoic_acid + 2H2O2 -> Azelaic_acid + Nonanoic_acid+ 2H2O', 'Nonanal', X = X2),
        #            ])
          # selectivity_Azelaic_acid = selectivity_Nonanoic_acid = X2 * X1
        # selectivity_Nonanal = selectivity_Oxonanoic_acid = X1 * (1 - X2)
        # selectivity_oxiraneoctanoic_acid,_3-octyl- = (1 - X1)
        
        #Considering that 
        # X=self.Oleic_acid_conversion
        # X1 = 1 - (self.selectivity_oxiraneoctanoic_acid)
        # X2 = self.selectivity_Azelaic_acid / X1

# class Separator(bst.Unit):
  
#     _N_outs = 6
        
#     def _run(self):
#         feed = self.ins[0]
#         IDs = ['Oleic_acid','Nonanal','Nonanoic_acid','Azelaic_acid', 
#                   'Oxononanoic_acid', 'oxiraneoctanoic_acid,_3-octyl-','Water']
#         outs = self.outs[0]
        
#         for stream, ID in zip(self.outs,IDs):
#             stream.imol[ID] = feed.imol[ID]
           












        
        
        
        
        
        
        
