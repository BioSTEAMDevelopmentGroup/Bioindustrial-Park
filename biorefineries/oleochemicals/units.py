"""
Created on Fri Oct 29 08:17:38 2021
@author: yrc2
"""
import biosteam as bst
import thermosteam as tmo
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream

class OxidativeCleavageReactor(bst.BatchBioreactor):
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
            tmo.Rxn('Oleic_acid + H2O2   -> Epoxy_stearic_acid + Water ', 'Oleic_acid', X=self.Oleic_acid_conversion),
            tmo.Rxn('Epoxy_stearic_acid + H2O2 -> Nonanal + Oxononanoic_acid + H2O', 'oxiraneoctanoic_acid,_3-octyl-', X = X1),
            tmo.Rxn('Nonanal + Oxononanoic_acid + 2H2O2 -> Azelaic_acid + Nonanoic_acid+ 2H2O', 'Nonanal', X = X2),
                    ])

        
        #NEEDS TO BE CONFIRMED WITH DANIM
        
        # Epoxide_formation = SRxn([Rxn('Oleic_acid + H2O2   -> Epoxy_stearic_acid ', 'Oleic_acid', X=self.Oleic_acid_conversion),
        #                           Rxn('Epoxy_stearic_acid + H2O -> DHSA', 'Epoxy_stearic_acid', X = 1),
        #                           Rxn('DHSA -> Nonanal + Oxononanoic_acid', 'DHSA', X = X1)])
        
        # Product_formation = PRxn([ Rxn('Nonanal + H2O2 ->  Nonanoic_acid', 'Nonanal', X = X2),
        #                            Rxn('Oxononanoic_acid + H2O2 -> Azelaic_acid', 'Oxononanoic_acid', X = X2),
        #                         ])

        # oxidative_cleavage_rxnsys = RxnSys(Epoxide_formation,Product_formation)
        # self.reactions = oxidative_cleavage_rxnsys
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)
              
        self.reactions(effluent) 
        
        
        effluent.T = self.T
        effluent.P = self.P
        
class Separator(bst.Unit):
  
    _N_outs = 6
        
    def _run(self):
        feed = self.ins[0]
        IDs = ['Oleic_acid','Nonanal','Nonanoic_acid','Azelaic_acid', 
                  'Oxononanoic_acid', 'oxiraneoctanoic_acid,_3-octyl-','Water']
        outs = self.outs[0]
        
        for stream, ID in zip(self.outs,IDs):
            stream.imol[ID] = feed.imol[ID]
           
            
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

        
        # self.reactions = tmo.SeriesReaction([
        #     tmo.Rxn('Oleic_acid + H2O2   -> Epoxy_stearic_acid + Water ', 'Oleic_acid', X=self.Oleic_acid_conversion),
        #     tmo.Rxn('Epoxy_stearic_acid + H2O2 -> Nonanal + Oxononanoic_acid + H2O', 'oxiraneoctanoic_acid,_3-octyl-', X = X1),
        #     tmo.Rxn('Nonanal + Oxononanoic_acid + 2H2O2 -> Azelaic_acid + Nonanoic_acid+ 2H2O', 'Nonanal', X = X2),
        #            ])
  
    
# #Solubility data
# #Pelargonic acid is insoluble in water
# #Interpolation for Azelaic acid in water
# # m = 22-2/50-20 = 1.96
# # S = mT + C
# # S = 1.96T - 76.0
# #Therefore, we assume azelaic acid solubility at 95 deg cel is 110g/l
















        
        
        
        
        
        
        
