# # -*- coding: utf-8 -*-
# """
# Created on Fri Oct 29 08:17:38 2021
# @author: yrc2
# """


# import biosteam as bst
# import thermosteam as tmo
 
#  class OzonolysisReactor(bst.BatchBioreactor):
#      _N_ins = 1
#      _N_outs = 1
     
#      @property
#      def effluent(self):
#          return self.outs[0]
     
#      @effluent.setter
#      def effluent(self,effluent):
#          self.outs[0]=effluent
         
     
#      def __init__(self, ID='', ins=None, outs=(), thermo=None,
#                   tau=17, N=None, V=None, T=373.15, P=101325,
#                   Nmin=2, Nmax=36):
         
#          bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
#                                     tau = tau , N = N, V = V, T = T, 
#                                     P = P ,Nmin = Nmin , Nmax = Nmax)
         
#  #When you inherit, you write all the parameters of the base class. Why is Nmin = Nmin
 
#          c = self.chemicals      
         
#          self.Oleic_acid_conversion = 0.808
#          #Oleic_acid formula = C18H34O2
#          #EPS formula = C18H34O3
         
#          #self.nonanal_selectivity = 0.094 
#          #mol EPS/ mol OA
                              
#          self.selectivity_oxiraneoctanoic_acid = 0.071
#          #self.Oxononanoic_acid_selectivity = 0.029
#          #'C9H16O3'
         
#          #self.Nonanoic_acid_selectivity = 0.277 
#          #C9H18O2
         
         
#          self.selectivity_Azelaic_acid = 0.442 * 2 
#          #carbon moles of EPS*(1 mol of EPS/ 9 C moles)/ carbon moles of oleic acid*(1 mol of OA/ 18 C moles)
#          # AA C9H16O4
         
#      def _setup(self):
#          # selectivity_Azelaic_acid = selectivity_Nonanoic_acid = X2 * X1
#          # selectivity_Nonanal = selectivity_Oxonanoic_acid = X1 * (1 - X2)
#          # selectivity_oxiraneoctanoic_acid,_3-octyl- = (1 - X1)
         
#          #Considering that 
#          X=self.Oleic_acid_conversion
#          X1 = 1 - (self.selectivity_oxiraneoctanoic_acid)
#          X2 = self.selectivity_Azelaic_acid / X1
         
#          self.reactions = tmo.SeriesReaction([
#              tmo.Rxn('Oleic_acid + H2O2 -> Epoxy_stearic_acid + Water ', 'Oleic_acid', X=self.Oleic_acid_conversion),
#              tmo.Rxn('Epoxy_stearic_acid + H2O2 -> Nonanal + Oxononanoic_acid + H2O', 'oxiraneoctanoic_acid,_3-octyl-', X = X1),
#              tmo.Rxn('Nonanal + Oxononanoic_acid + 2H2O2 -> Azelaic_acid + Nonanoic_acid+ 2H2O', 'Nonanal', X = X2)
#              ])
         
#      def _run(self):
#          feed = self.ins[0]
#          effluent = self.outs[0]
         
#          #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
#          effluent.copy_like(feed)
               
#          self.reactions(effluent) 
         
#          effluent.T = self.T
#          effluent.P = self.P
         
 
 
#  class Separator(bst.Unit):
#      #Why does this not have an _init_?
     
#      _N_outs = 6
         
#      def _run(self):
#          feed = self.ins[0]
#          IDs = ['Oleic_acid','Nonanal','Nonanoic_acid','Azelaic_acid', 
#                    'Oxononanoic_acid', 'oxiraneoctanoic_acid,_3-octyl-','Water']
#          outs = self.outs[0]
         
#          for stream, ID in zip(self.outs,IDs):
#              stream.imol[ID] = feed.imol[ID]
             
 
#  class AACrystalliser(bst.BatchCrystallizer):
     
#      def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
#                   T, crystal_AA_purity=0.95
#                   ):
#          bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
#                                         tau=2, V=1e6, T=T)
#          self.crystal_AA_purity = crystal_AA_purity
         
         
#      @property
#      def Hnet(self):
#          feed = self.ins[0]
#          effluent = self.outs[0]
#          if 's' in feed.phases:
#              H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, feed['s'].mol) if i.Hfus])
#          else:
#              H_in = 0.
#          solids = effluent['s']
#          H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
#          return H_out - H_in
         
#      def _run(self):
#          outlet = self.outs[0]
#          outlet.phases = ('s', 'l')
#          crystal_AA_purity = self.crystal_AA_purity
#          feed = self.ins[0]
#          AA_solid = crystal_AA_purity * outlet.imass['s']
#          rest_solid = (outlet.imass['s'] - AA_solid)/(1 - crystal_AA_purity)
#          AA_liquid = feed.imass['Azelaic_acid'].value - AA_solid
#          rest_liquid = outlet.imass['l'] - AA_liquid
         
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
        
class Separator(bst.Unit):
    #Why does this not have an _init_?
    
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
                        
    @property
    def Hnet(self):
        effluent = self.outs[0]
        solids = effluent['s']
        c = solids.chemicals
        H_tot = - (
            c.Azelaic_acid.Hf * AA_solid + c.Nonanoic_acid.Hf * NA
        )                  
           
        return H_tot
        
    def solubility(self, T):        
        m = (self.AA_molefraction_330_15K - self.AA_molefraction_280_15K) / (330.15 - 280.15)
        c = (self.AA_molefraction_330_15K - (m* 330.15))
        S = (m*(self.T) + c)
        return S
         

#Assuming inlet at saturation.
    def _run(self):
        outlet = self.outs[0]
        outlet.phases = ('s', 'l')
        feed = self.ins[0]    
        outlet.copy_like(feed)
        outlet.T = self.T
        x = self.solubility(outlet.T)
        outlet.sle('Azelaic_acid', solubility=x)
        # AA_solid = outlet.imol['s',('Azelaic_acid')]
        # NA = feed.imol ['Nonanoic_acid']
        
        
                

#         feed = self.ins[0]
#         AA, NA, Hexane, Water = feed.imass['Azelaic_acid',
#                                            'Nonanoic_acid',
#                                            'Hexane',
#                                            'Water'].value
#         total = AA + NA  
#         outlet.empty()
#         AA_solid_molefrac = 1 - self.AA_molefraction_280K         
#         AA_in_liquid = (self.AA_molefraction_280K *feed.imass['Water'].value * 188.8)/( AA_solid_molefrac * 18)
# # NA_in_liquid = self.NA_solubility*feed.imass['Water'] if only cooling till 20 deg cel
# # Assumption is that all the NA crystallises at 280K
#         AA_solid = AA - AA_in_liquid
#         outlet.imass['s',('Azelaic_acid','Nonanoic_acid')] =  [AA_solid, NA] 
#         outlet.imass['l',('Azelaic_acid','Hexane','Water')] = [AA_in_liquid,Hexane,Water]
#         outlet.T = self.T
      

        
        
        
        
        
        
        

