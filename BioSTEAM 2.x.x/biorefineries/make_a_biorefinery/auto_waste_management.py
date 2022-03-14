# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 19:01:59 2021

@author: sarangbhagwat
"""

from biosteam import Facility

# from . import Facility

__all__ = ('AutoWasteManagement')


class AutoWasteManagement(Facility):
    """
    Create an AutoWasteManagement object. 
    When run, all streams with IDs containing the String to_boiler_solids_mixer_ID_key are mixed as the final inlet of wastewater_mixer, and
    all streams with IDs containing the String to_boiler_solids_mixer_ID_key are mixed as the final inlet of boiler_solids_mixer.
    If any streams with those labels were manually set as an inlet for either mixer, those streams will be excluded from the final inlet 
    of that mixer. Each mixer must be defined with its final (or only) having an ID string.
    
    Parameters
    ----------
    wastewater_mixer: Unit
        Used by the biorefinery to mix wastes to be sent to wastewater treatment and/or disposal.
    
    boiler_solids_mixer: Unit
        Used by the biorefinery to mix wastes to be sent to the boiler for combustion.
    
    to_wastewater_mixer_ID_key: str, optional
        Used to identify streams in the flowsheet to be sent to the wastewater_mixer unit.
        Defaults to '_to_WWT'.
        
    to_boiler_solids_mixer_ID_key: str, optional
        Used to identify streams in the flowsheet to be sent to the boilder_solids_mixer unit.
        Defaults to '_to_boiler'
        
    Examples
    -------
    >>>import thermosteam as tmo
    >>>import biosteam as bst
    >>>from biosteam import SystemFactory
    >>>flowsheet = bst.Flowsheet('AWM_test')
    >>>bst.main_flowsheet.set_flowsheet(flowsheet)
    
    >>>Lignin=tmo.Chemical('Lignin')
    >>>Lignin.Psat.add = 0.
    >>>Lignin.Tb = 600
    >>>Lignin.Hvap = 1e5
    >>>tmo.settings.set_thermo([Lignin, 'Water', 'AceticAcid'])
    
    >>>@SystemFactory(ID = 'AWM_test_sys')
    >>>def create_AWM_test_sys(ins, outs):
    >>>    S101_feed=tmo.Stream('S101_feed', Lignin=10, Water=30, AceticAcid=10)
    >>>    M101_feed=tmo.Stream('M101_feed', Lignin=5)
    >>>    S101 = bst.units.FakeSplitter('S101', ins=S101_feed, outs=('S101_to_WWT', 'S101_to_boiler'))
    >>>    @S101.add_specification()
    >>>    def S101_spec():
    >>>        S101.outs[1].imol['Lignin'] = S101.ins[0].imol['Lignin']
    >>>        S101.outs[0].copy_like(S101.ins[0])
    >>>        S101.outs[0].imol['Lignin'] -= S101.ins[0].imol['Lignin']
    >>>    M101 = bst.units.Mixer('M101', ins=M101_feed, outs=('M101_to_boiler'))
    >>>    M501 = bst.units.Mixer('M501', ins=('mixed_wastewater_to_M501'), outs=())
    >>>    M505 = bst.units.Mixer('M505', ins=(M101-0,'mixed_solids_to_M505'), outs=())
    >>>    AWM = bst.facilities.AutoWasteManagement('AWM', wastewater_mixer=M501, boiler_solids_mixer=M505,
    >>>                              to_wastewater_mixer_ID_key='to_WWT',
    >>>                              to_boiler_solids_mixer_ID_key='to_boiler')
    >>>AWM_test_sys = create_AWM_test_sys()
    >>>AWM_test_sys.simulate()
    >>>u = flowsheet.unit
    >>>u.M501.show()
    >>>u.M505.show()
    

    The AutoWasteManagement facility is being used. Note that
    all streams labeled 'to_WWT' are mixed as the final inlet of M501, and
    all streams labeled 'to_boiler' are mixed as the final inlet of M505.
    If any streams with those labels were manually set as an inlet for either mixer, those
    streams will be excluded from the mixed final inlet of that mixer.
    
    
    Mixer: M501
    ins...
    [0] mixed_wastewater_to_M501
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water       30
                        AceticAcid  10
    outs...
    [0] s10
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water       30
                        AceticAcid  10
    Mixer: M505
    ins...
    [0] M101_to_boiler  from  Mixer-M101
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Lignin  5
    [1] mixed_solids_to_M505
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Lignin  10
    outs...
    [0] s12
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Lignin  15
    """
    network_priority = 100 # must be prior to wastewater treatment and boiler
    _N_ins = 0
    _N_outs = 0
    _N_heat_utilities = 0
    _N_power_utilities = 0
    
    def __init__(self, ID='',
                 wastewater_mixer=None, boiler_solids_mixer=None,
                 to_wastewater_mixer_ID_key='_to_WWT',
                 to_boiler_solids_mixer_ID_key='_to_boiler'):
        Facility.__init__(self, ID, ins=None, outs=None, thermo=None)
        self.wastewater_mixer = wastewater_mixer
        self.boiler_solids_mixer = boiler_solids_mixer
        self.to_wastewater_mixer_ID_key = to_wastewater_mixer_ID_key
        self.to_boiler_solids_mixer_ID_key = to_boiler_solids_mixer_ID_key
        
        print("\n\nThe AutoWasteManagement facility is being used. Note that")
        print(f"all streams labeled '{to_wastewater_mixer_ID_key}' are mixed as the final inlet of {wastewater_mixer.ID}, and")
        print(f"all streams labeled '{to_boiler_solids_mixer_ID_key}' are mixed as the final inlet of {boiler_solids_mixer.ID}.")
        print("If any streams with those labels were manually set as an inlet for either mixer, those")
        print("streams will be excluded from the mixed final inlet of that mixer (i.e., the manual connection will be retained.\n\n")
        
        @wastewater_mixer.add_specification(run=True)
        def wastewater_mixer_spec():
            self._run()
            
        @boiler_solids_mixer.add_specification(run=True)
        def boiler_solids_mixer_spec():
            self._run()
                
    def _run(self):
        wastewater_mixer = self.wastewater_mixer
        boiler_solids_mixer = self.boiler_solids_mixer
        to_wastewater_mixer_ID_key = self.to_wastewater_mixer_ID_key
        to_boiler_solids_mixer_ID_key = self.to_boiler_solids_mixer_ID_key

        s = self.system.flowsheet.stream
        
        self.individual_streams_to_WWT = individual_streams_to_WWT = []
        self.individual_solids_to_boiler = individual_solids_to_boiler = []
        
        self.mixed_streams_to_WWT = mixed_streams_to_WWT = wastewater_mixer.ins[-1]
        self.mixed_solids_to_boiler = mixed_solids_to_boiler = boiler_solids_mixer.ins[-1]
        
        for si in s:
            if to_wastewater_mixer_ID_key in si.ID and \
                not (si in wastewater_mixer.ins or si in wastewater_mixer.outs) and \
                not si.source in wastewater_mixer.get_downstream_units():
                individual_streams_to_WWT.append(si)
            elif to_boiler_solids_mixer_ID_key in si.ID and \
                not (si in boiler_solids_mixer.ins or si in boiler_solids_mixer.outs) and \
                not si.source in boiler_solids_mixer.get_downstream_units():
                individual_solids_to_boiler.append(si)
        
        mixed_streams_to_WWT.empty()
        mixed_streams_to_WWT.mix_from(individual_streams_to_WWT)
        mixed_solids_to_boiler.empty()
        mixed_solids_to_boiler.mix_from(individual_solids_to_boiler)
        
        self._run_wastewater_mixer()
        self._run_boiler_solids_mixer()
        
    def _design(self): pass
    
    def _cost(self): pass

    @property
    def wastewater_streams(self):
        return self.wastewater_mixer.ins
    @property
    def solid_waste_streams(self):
        return self.boiler_solids_mixer.ins
    
    @property
    def mixed_wastewater(self):
        return self.wastewater_mixer.outs[0]
    @property
    def mixed_solid_waste(self):
        return self.boiler_solids_mixer.outs[0]
    
    def _run_wastewater_mixer(self):
        self._cost()
        self.wastewater_mixer._run()
    
    def _run_boiler_solids_mixer(self):
        self._cost()
        self.boiler_solids_mixer._run()
        