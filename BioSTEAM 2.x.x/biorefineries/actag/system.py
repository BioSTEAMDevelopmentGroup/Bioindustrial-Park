# -*- coding: utf-8 -*-
"""
The complete acTAG biorefinery system is created here.

"""
import flexsolve as flx
import thermosteam as tmo
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import SystemFactory
from ..cornstover import (
    create_ammonia_fiber_expansion_pretreatment_system,
    create_cellulosic_fermentation_system,
)
from ..lipidcane import (
    create_feedstock_handling_system,
    create_juicing_system,
    create_lipid_pretreatment_system as create_oil_pretreatment_system,
    create_transesterification_and_biodiesel_separation_system,
)
import biorefineries as brf
from biorefineries.oilcane import units
from thermosteam import Rxn, PRxn

__all__ = (
    'create_cellulosic_acTAG_system',
    'create_conventional_acTAG_system',
)

@bst.SystemFactory(
    ID='cellulosic_acTAG_sys',
    ins=[dict(ID='switchgrass', 
              Arabinan=0.02789023841655421,
              Galactan=0.010436347278452543,
              Glucan=0.2717049032838507,
              Xylan=0.21214574898785432,
              Mannan=0.005937921727395412,
              Lignin=0.17112010796221322,
              Ash=0.016194331983805668,
              Extractives=0.08457040035987407,
              Water=0.2,
              total_flow=104229.16,
              units='kg/hr')],
    outs=[dict(ID='acTAG')],
)
def create_cellulosic_acTAG_system(ins, outs):
    feedstock, = ins
    acTAG, = outs
    chemicals = feedstock.chemicals
    AFEX_sys = create_ammonia_fiber_expansion_pretreatment_system(
        feedstock,
        include_feedstock_handling=True,
        solids_loading=0.625,
        ammonia_loading=0.555,
        T_pretreatment_reactor=273.15 + 100.,
        residence_time=0.5,
        pretreatment_reactions=PRxn([
        #            Reaction definition                 Reactant    Conversion
        Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.0050, chemicals),
        Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030, chemicals),
        Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.0030, chemicals),
        Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1.0000, chemicals),
        Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9000, chemicals),
        Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.0500, chemicals),
        Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9000, chemicals),
        Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0050, chemicals),
        Rxn('Acetate -> AceticAcid',                     'Acetate',  1.0000, chemicals),
        Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0500, chemicals)
            ]),
    )
    hydrolysate, pretreatment_wastewater = AFEX_sys.outs
    
    

@SystemFactory(
    ID='conventional_acTAG_sys',
    ins=[dict(ID='sugarcane',
              Water=0.7,
              Glucose=0.01208,
              Sucrose=0.1369,
              Ash=0.006,
              Cellulose=0.06115,
              Hemicellulose=0.03608,
              Lignin=0.03276,
              Solids=0.015,
              total_flow=333334.2,
              units='kg/hr',
              price=0.03455)],
    outs=[dict(ID='acTAG')]
)
def create_conventional_acTAG_system(ins, outs):
    feedstock, = ins
    acTAG, = outs
    create_feedstock_handling_system,
    create_juicing_system,