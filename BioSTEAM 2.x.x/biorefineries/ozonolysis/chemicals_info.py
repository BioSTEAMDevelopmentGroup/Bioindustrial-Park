# -*- coding: utf-8 -*-
"""
"""
import thermosteam as tmo
import biosteam as bst
from biosteam import Unit, Stream, settings, main_flowsheet
#%% Chemicals - Definition
# Create chemicals here
ozo_chemicals = tmo.Chemicals(
    ['Water','Hydrogen_peroxide','Oleic_acid',
     'Nonanal','Nonanoic_acid','Azelaic_acid'])

(Water,Hydrogen_peroxide,Oleic_acid,
Nonanal,Nonanoic_acid,Azelaic_acid) = ozo_chemicals


def create_new_chemical(ID, phase='s', **constants):
    # Create a new solid chemical without any data
    solid = tmo.Chemical(ID, search_db=False, phase=phase, **constants)

    # Add chemical to the Chemicals object
    ozo_chemicals.append(solid)

    return solid
#Writing chemicals not defined in the database and other additional chemicals

Catalyst = create_new_chemical(
    'Phosphotungstic_acid',
    formula="H3PW12O40", # Chemical Formula
    MW=2880.2,
    CAS='1343-93-7'
    )

Oxononanoic_acid = create_new_chemical('Oxononanoic_acid',
    phase='l',
    Hf=-579480,
    formula = 'C9H16O3',
    MW = 172.22,
    CAS = '2553-17-5'
)

Epoxide = create_new_chemical('Epoxy_stearic_acid',
    phase='l',
    #Hf,
    formula = 'C18H34O3',
    MW = 298.5,
    CAS = '2443-39-2'
 )
for chemical in ozo_chemicals: chemical.default()

tmo.settings.set_thermo(ozo_chemicals)
 




#%% Stream Data and Mass balance
mixed_feed_stream = tmo.Stream('mixed_feed_stream')
mixed_feed_stream.imol['Oleic_acid']=0.86
mixed_feed_stream.imol['H2O2']=6.85
mixed_feed_stream.imol['H2O']=27.1
print(mixed_feed_stream.F_mass)
# Ozonolysis_series_rxn(mixed_feed_stream)
# print(mixed_feed_stream.F_mass)
# mixed_feed_stream.show(N=100)

# outs = [Stream('reactor_out')]
#!!!TODO
#Change conversion values
#Check if you need diol data

#%% Units

import biosteam as bst

class OzonolysisReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
    @property
    def effluent(self):
        return self.outs[0]

    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=17, N=None, V=None, T=373.15, P=101325,
                 Nmin=2, Nmax=36):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                   tau = tau , N = N, V = V, T = T, 
                                   P =P ,Nmin =Nmin, Nmax = Nmax)
    
        
        
    def _setup(self):
        self.reactions = Ozonolysis_parallel_rxn = tmo.SeriesReaction([
             #Assumption, every conversion is 0.947 and overall conversion is (0.947^3)
            tmo.Rxn('Oleic_acid + H2O2 -> Epoxy_stearic_acid+ Water ', 'Oleic_acid', X= 0.947),
            tmo.Rxn('Epoxy_stearic_acid + H2O2 -> Nonanal + Oxononanoic_acid + H2O', 'Epoxy_stearic_acid', X = 0.947),
            tmo.Rxn('Nonanal + Oxononanoic_acid + 2H2O2 -> Azelaic_acid + Nonanoic_acid+ 2H2O', 'Nonanal', X = 0.947)]
        )
        #Ozonolysis_parallel_rxn.correct_atomic_balance(['Oleic_acid','H2O2','9_10_epoxy_stearic_acid','H2O','Nonanal','9_Oxononanoic_acid','Nonaoic_acid','Azelaic_acid'])
    
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        effluent.copy_like(feed)
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        
reactor = OzonolysisReactor(ID ='',ins = mixed_feed_stream)
reactor.simulate()
print(reactor.results())
reactor.show()

 
     





