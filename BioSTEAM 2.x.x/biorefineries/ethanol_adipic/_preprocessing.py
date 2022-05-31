#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

import pandas as pd
import biosteam as bst
from biosteam import Unit
from biosteam.utils.misc import format_title
from biosteam.units.decorators import add_cost
from biorefineries.lactic._utils import CEPCI, dry_composition

__all__ = ('prep_cost', 'Preprocessing', 'OptionPrep',
           'Auxiliary', 'Grinder', 'HammerMill', 'Dryer', 'PelletMill', 'DepotAFEX',
           'create_default_depot', 'PreprocessingCost')


# %%

def prep_cost(basis, ID=None, *, CE, cost, n,
              lifetime, salvage, maintenance, labor,
              S=1, lb=None, ub=None, kW=0, BM=1, units=None, f=None, N=None):

    def add_param(cls):
        add_cost(cls, ID, basis, units, S, lb, ub, CE, cost, n, kW, BM, N, lifetime, f)
        cls.salvage[ID] = salvage
        cls.maintenance[ID] = maintenance
        cls.labor[ID] = labor
        return cls

    return add_param


class Preprocessing(Unit):
    '''
    For preprocessing units at biomass depots based on description in Lamers et al., [1]_
    capacity of the plant is set at 9.07 metric tonne (MT, or Mg) per hour
    (equivalent to 10 U.S. ton per hr).

    Parameters
    ----------
    ins : biosteam.Stream
        [0] Feedstock
    outs : biosteam.Stream
        [0] Feedstock
    salvage : dict
        Salvage value as a fraction of purchase cost for each equipment in this unit.
    maintenance : dict
        Maintenance cost as a fraction of purchase cost for each equipment in this unit.
    labor : dict
        Labor cost per hour for each equipment in this unit.

    Note
    ----
    It is more convenient to use the prep_cost decorator to add equipment information.

    References
    ----------
    .. [1] Lamers et al., Techno-Economic Analysis of Decentralized Biomass
        Processing Depots. Bioresource Technology 2015, 194, 205–213.
        https://doi.org/10.1016/j.biortech.2015.07.009.

    '''

    salvage = {}
    maintenance = {}
    labor = {}

    _units= {'Dry flow': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), **kwargs):
        Unit.__init__(self, ID, ins, outs)
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def _run(self):
        self.outs[0].copy_like(self.ins[0])

    def _design(self):
        self.design_results['Dry flow'] = self.ins[0].F_mass - self.ins[0].imass['Water']

    @property
    def moisture_in(self):
        '''[float] Moisture content of feedstock prior to this unit'''
        return self.ins[0].imass['Water']/self.ins[0].F_mass

    @property
    def moisture_out(self):
        '''[float] Moisture content of feedstock after this unit'''
        return self.outs[0].imass['Water']/self.outs[0].F_mass


class OptionPrep(Preprocessing):

    def __init__(self, ID='', ins=None, outs=(), kind='CPP', **kwargs):
        Unit.__init__(self, ID, ins, outs)
        self.kind = kind
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def _cost(self):
        name = format_title(type(self).__name__)
        ID = f'{self.kind} {name.lower()}'
        self.cost_items = {ID: self.cost_items[ID]}
        self.purchase_costs.clear()
        self._decorated_cost()

    @property
    def kind(self):
        '''
        [str] Type of the equipment, can be either CPP or HMPP (conventional or
        high-moisture pelleting process).
        '''
        return self._kind
    @kind.setter
    def kind(self, i):
        if i not in ('CPP', 'HMPP'):
            raise ValueError(f'kind can only be CPP or HMPP, not {i}.')
        self._kind = i



# %%


@prep_cost(basis='Dry flow', ID='Conveyor', units='kg/hr',
           kW=0.167*9.07, cost=268800, S=9070, CE=CEPCI[2011], n=1, BM=1,
           lifetime=168000/(365*24), salvage=0.3, maintenance=0,
           labor=0)
@prep_cost(basis='Dry flow', ID='Dust collection', units='kg/hr',
           kW=9.33*9.07, cost=286400, S=9070, CE=CEPCI[2011], n=1, BM=1,
           lifetime=168000/(365*24), salvage=0.3, maintenance=0,
           labor=0)
@prep_cost(basis='Dry flow', ID='Surge bin', units='kg/hr',
           kW=0.167*9.07, cost=96800, S=9070, CE=CEPCI[2011], n=1, BM=1,
           lifetime=168000/(365*24), salvage=0.3, maintenance=0,
           labor=0)
@prep_cost(basis='Dry flow', ID='Miscellaneous', units='kg/hr',
           kW=0.333*9.07, cost=84000, S=9070, CE=CEPCI[2011], n=1, BM=1,
           lifetime=168000/(365*24), salvage=0.3, maintenance=0,
           labor=0)
class Auxiliary(Preprocessing):
    '''
    Auxiliary units including conveyor system, dust collection system, surge bin,
    and other miscellaneous equipment.

    Note
    ----
    [1] For conveyor, dust collection, surge bin, and miscellaneous equipment,
        lifetime, salvage, and electricity usage was not explicitly provided in
        Lamers et al., [1]_ thus was back calculated based on the provided information.

    References
    ----------
    .. [1] Lamers et al., Techno-Economic Analysis of Decentralized Biomass
        Processing Depots. Bioresource Technology 2015, 194, 205–213.
        https://doi.org/10.1016/j.biortech.2015.07.009.

    '''

@prep_cost(basis='Dry flow', ID='Grinder', units='kg/hr',
           kW=19.85*4.54, cost=180000, S=4540, CE=CEPCI[2011], n=1, BM=1,
           lifetime=15000/(365*24), salvage=0.3, maintenance=0.1,
           labor=0.5*15.88)
class Grinder(Preprocessing): pass


@prep_cost(basis='Dry flow', ID='CPP hammer mill', units='kg/hr',
           kW=13.23*4.54, cost=105225, S=4540, CE=CEPCI[2011], n=1, BM=1,
           lifetime=40000/(365*24), salvage=0.3, maintenance=0.1,
           labor=0.5*19.88)
@prep_cost(basis='Dry flow', ID='HMPP hammer mill', units='kg/hr',
           kW=66.15*1.81, cost=105225, S=1810, CE=CEPCI[2011], n=1, BM=1,
           lifetime=40000/(365*24), salvage=0.3, maintenance=0.1,
           labor=0.5*19.88)
class HammerMill(OptionPrep):

    _cost = OptionPrep._cost


@prep_cost(basis='Dry flow', ID='CPP dryer', units='kg/hr',
           kW=330.76*1.81, cost=354310, S=1810, CE=CEPCI[2011], n=1, BM=1,
           lifetime=168000/(365*24), salvage=0.3, maintenance=0.1,
           labor=0.5*15.51)
@prep_cost(basis='Dry flow', ID='HMPP dryer', units='kg/hr',
           kW=55.13*4.54, cost=35009, S=4540, CE=CEPCI[2011], n=1, BM=1,
           lifetime=168000/(365*24), salvage=0.3, maintenance=0.1,
           labor=0.5*15.51)
class Dryer(OptionPrep):

    _target_moisture = 0.2

    def _run(self):
        mc_out = self.target_moisture
        dried = self.outs[0]
        dried.copy_like(self.ins[0])
        dry_mass = dried.F_mass - dried.imass['Water']
        dried.imass['Water'] = mc_out/(1-mc_out) * dry_mass

    def _cost(self):
        OptionPrep._cost(self)
        mc_diff = self.moisture_in - self.moisture_out
        if mc_diff < 0:
            self.power_utility.rate = 0
            self.purchase_costs = dict.fromkeys(self.purchase_costs.keys(), 0)
        elif self.kind == 'CPP':
            self.power_utility.rate *= mc_diff/(0.3-0.09)
        else:
            self.power_utility.rate *= mc_diff/(0.3-0.19)

    @property
    def target_moisture(self):
        '''[float] Target moisture content after drying.'''
        return self._target_moisture
    @target_moisture.setter
    def target_moisture(self, i):
        self._target_moisture = float(i)


@prep_cost(basis='Dry flow', ID='CPP pellet mill', units='kg/hr',
           kW=55.13*4.54, cost=304264, S=4540, CE=CEPCI[2011], n=1, BM=1,
           lifetime=100000/(365*24), salvage=0.3, maintenance=0.1,
           labor=0.5*15.51)
@prep_cost(basis='Dry flow', ID='HMPP pellet mill', units='kg/hr',
           kW=82.69*4.54, cost=304264, S=4540, CE=CEPCI[2011], n=1, BM=1,
           lifetime=100000/(365*24), salvage=0.3, maintenance=0.1,
           labor=0.5*15.51)
class PelletMill(OptionPrep):

    def _run(self):
        if self.kind == 'CPP':
            OptionPrep._run(self)
        else:
            mc_out = self.moisture_in - 0.1
            dried = self.outs[0]
            dried.copy_like(self.ins[0])
            dry_mass = dried.F_mass - dried.imass['Water']
            dried.imass['Water'] = mc_out/(1-mc_out) * dry_mass

    _cost = OptionPrep._cost


# J/mol / (g/mol) = J/g = kJ/kg
CH4_HHV = bst.Chemical('CH4').HHV/bst.Chemical('CH4').MW
@prep_cost(basis='Dry flow', ID='AFEX', units='kg/hr',
           kW=60.64*8.64, cost=3101027, S=8640, CE=CEPCI[2011], n=1, BM=1,
           lifetime=262800/(365*24), salvage=0.15, maintenance=0.05,
           labor=1*19.88)
class DepotAFEX(Preprocessing):
    '''
    Ammonia fiber expansion (AFEX) pretreatment at a biorefinery depot,
    equipment information based on description in Lamers et al., [1]_
    makeup water and ammonia (makeup and fugative) based on Kim and Dale. [2]_

    Parameters
    ----------
    ins : Stream sequence
        [0] Feedstock
        [1] Makeup water
        [2] Makeup ammonia for AFEX pretreatment, note that this stream should have price,
            but its flow will be updated based on the given loading
        [3] Natural gas for heating, note that this stream should have price,
            but its flow will be updated based on the given loading
    outs : Stream sequence
        [0] Feedstock
        [1] Fugative ammonia
    loss : float
        Feedstock dry matter loss due to AFEX pretreatment.
    ammonia_makeup : float
        Makeup ammonia as kg ammonia/kg dry feedstock.
    ammonia_fugative : float
        Fugative ammonia as kg ammonia/kg dry feedstock.
    CH4_loading : float
        Natural gas consumption as kg CH4/kg dry feedstock.

    References
    ----------
    .. [1] Lamers et al., Techno-Economic Analysis of Decentralized Biomass
        Processing Depots. Bioresource Technology 2015, 194, 205–213.
        https://doi.org/10.1016/j.biortech.2015.07.009.
    .. [2] Kim and Dale., Comparing Alternative Cellulosic Biomass Biorefining
        Systems: Centralized versus Distributed Processing Systems.
        Biomass and Bioenergy 2015, 74, 135–147.
        https://doi.org/10.1016/j.biombioe.2015.01.018.

    '''

    _N_ins = 4
    _N_outs = 2

    loss = 0
    water_makeup = 2.1/3.6
    ammonia_makeup = 31.3/1e3/3.6
    ammonia_fugative = 11.1/1e3/3.6
    # 325 kWh/Mg as in ref [1], did not use value in ref [2] as ref [2] also has
    # steam consumption
    CH4_loading = (325*3600/-CH4_HHV)/1e3

    def _run(self):
        feed, makeup_H2O, makeup_NH3, CH4  = self.ins
        processed, fugative_NH3 = self.outs

        dry_mass = feed.F_mass - feed.imass['Water']
        makeup_H2O.imass['Water'] = self.water_makeup * dry_mass
        makeup_NH3.imass['NH3'] = self.ammonia_makeup * dry_mass
        CH4.imass['CH4'] = self.CH4_loading * dry_mass
        fugative_NH3.imass['NH3'] = self.ammonia_fugative * dry_mass
        fugative_NH3.phase = 'g'

        processed.copy_like(feed)
        processed.imol['NH4OH'] = makeup_NH3.imol['NH3'] - fugative_NH3.imol['NH3']
        processed.imol['H2O'] -= processed.imol['NH4OH']
        processed.F_mass *= 1 - self.loss


# %%

def create_default_depot(raw_feedstock='raw_feedstock', kind='CPP',
                         with_AFEX=False, water_price=0, ammonia_price=0,
                         natural_gas_price=0., prefix=1):
    '''
    Create a biorefinery depot system with default settings based on decriptions
    in Lamers et al., [1]_ note that the target moisture of the preprocessed
    feedstock is adjusted to 20% (as opposed to the 9% in [1]_) as described
    in later biorefinery reports. [2]_

    Parameters
    ----------
    raw_feedstock : biosteam.Stream or str
        Feedstock stream prior to preprocessing, if provided as a str, the stream
        will be automatically created based on the composition decribed in [2]_
        and the provided str is used as the ID.
    kind : str
        Either 'CPP' for conventional pelleting process or 'HMPP' for high-moisture
        pelleting process.
    with_AFEX : bool
        Whether to include ammonia fiber expansion pretreatment in the process.
    water_price : float
        Price of makeup water in AFEX, will only be used if with_AFEX=True.
    ammonia_price : float
        Price of makeup ammonia in AFEX, will only be used if with_AFEX=True.
    natural_gas_price : float
        Price of natural gas for heating in AFEX, will only be used if with_AFEX=True.
    prefix : int
        Prefix number for unit IDs.

    Returns
    -------
    flowsheet : biosteam.Flowsheet
        Flowsheet containing the depot system, the preprocessed feedstock can
        be retrieved as flowsheet.stream.preprocessed.

    References
    ----------
    .. [1] Lamers et al., Techno-Economic Analysis of Decentralized Biomass
        Processing Depots. Bioresource Technology 2015, 194, 205–213.
        https://doi.org/10.1016/j.biortech.2015.07.009.
    .. [2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
        Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
        NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
        https://doi.org/10.2172/1483234

    '''

    flowsheet = bst.Flowsheet(f'{kind}_depot')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    feed_ID = 'raw_feedstock'
    if isinstance(raw_feedstock, str):
        chems = bst.settings.get_chemicals()
        feed_ID = raw_feedstock
        dry_mass = 2000*1e3/24 # 2000 dry Mg per day
        mc = 0.3 # moisture content
        dry_array = chems.kwarray(dry_composition)
        wet_mass = dry_mass / (1-mc)
        moisture_array = chems.kwarray(dict(Water=mc))
        feedstock_flow = wet_mass * (dry_array*(1-mc)+moisture_array)
        raw_feedstock = bst.Stream(feed_ID, feedstock_flow, units='kg/hr')
    elif isinstance(raw_feedstock, bst.Stream):
        chems = raw_feedstock.chemicals
        bst.settings.set_thermo(chems)
    else:
        raise ValueError('raw_feedstock must be a Stream, str, or leave as blank, '
                          f'cannot be {type(raw_feedstock).__name__}.')

    n = int(prefix*100)
    U101 = Grinder(f'U{n+1}', ins=raw_feedstock)
    if with_AFEX:
        water = bst.Stream('water', price=water_price, units='kg/hr')
        ammonia = bst.Stream('ammonia', price=ammonia_price, units='kg/hr')
        natural_gas = bst.Stream('natural_gas', price=natural_gas_price,
                                 units='kg/hr')
        AFEX = DepotAFEX(f'U{n+2}', ins=(U101-0, water, ammonia, natural_gas))
        n += 1
    if kind == 'CPP':
        U102 = Dryer(f'U{n+2}', kind='CPP', target_moisture=0.2)
        U103 = HammerMill(f'U{n+3}', kind='CPP', ins=U102-0)
        U104 = PelletMill(f'U{n+4}', kind='CPP', ins=U103-0)
        U105 = Auxiliary(f'U{n+5}', ins=U104-0, outs='preprocessed')

    elif kind == 'HMPP':
        U102 = HammerMill(f'U{n+2}', kind='HMPP')
        U103 = PelletMill(f'U{n+3}', kind='HMPP', ins=U102-0)
        U104 = Dryer(f'U{n+4}', kind='HMPP', target_moisture=0.2, ins=U103-0)
        U105 = Auxiliary(f'U{n+5}', ins=U104-0, outs='preprocessed')

    else:
        raise ValueError(f'kind can only be CPP or HMPP, not {kind}.')

    if with_AFEX:
        AFEX-0-U102
        prep_sys = bst.System('prep_sys',
                              path=(U101, AFEX, U102, U103, U104, U105))
    else:
        U101-0-U102
        prep_sys = bst.System('prep_sys',
                              path=(U101, U102, U103, U104, U105))
    return flowsheet


# %%

class PreprocessingCost:
    '''
    To calculate the cost of feedstock preprocessing for a biorefinery depot
    ststen based on description in Lamers et al. [1]_

    Parameters
    ----------
    depot_sys : biosteam.System
        Depot system.
    feedstock : biosteam.Stream
        Feedstock stream of the depot system. Can leave as blank if the depot system
        only has one product stream.
    operating_hours : float
        Work hours per year.
    interest : float
        Annual interest rate.
    insurance : float
        Insurance cost as a fraction of purchase cost.
    housing : float
        Housing cost as a fraction of purchase cost.
    tax : float
        Taxes as a fraction of purchase cost.
    labor_adjustment : float
        Factor (relative to year 2011) for labor cost conversion.

    References
    ----------
    .. [1] Lamers et al., Techno-Economic Analysis of Decentralized Biomass
        Processing Depots. Bioresource Technology 2015, 194, 205–213.
        https://doi.org/10.1016/j.biortech.2015.07.009.

    '''

    __slots__ = ('ID', '_depot_sys', '_feedstock', '_units', '_throughput', '_feeds',
                 'operating_hours', 'interest', 'insurance',
                 'housing', 'tax', 'labor_adjustment')

    def __repr__(self):
        return f'<{type(self).__name__}: {self.depot_sys.ID}>'


    def __init__(self, depot_sys=None, feedstock=None,
                 operating_hours=365*24, interest=0.06,
                 insurance=0.0025, housing=0.0075, tax=0.01,
                 labor_adjustment=1):

        self.depot_sys = depot_sys
        if not feedstock:
            feedstock, = (i for i in depot_sys.products if i.phase != 'g')
        self.feedstock = feedstock
        self.operating_hours = 365*24
        self.interest = interest
        self.insurance = insurance
        self.housing = housing
        self.tax = tax
        self.labor_adjustment = labor_adjustment

    def get_depreciation_cost(self):
        r'''
        Calculate the depreciation cost per metric tonne (MT) of dry feedstock.

        .. math:: cost[\frac{$}{yr}] = (P-S) * \frac{i*(1+i)^n}{(1+i)^n-1} + S*i
        .. math:: cost[\frac{$}{MT}] = \frac{cost[\frac{$}{yr}]}{MT \ feedstock \ per \ yr}


        where:
            P: purchase cost

            S: salvage value = salvage fraction * purchase cost

            i: interest rate

            n: equipment lifetime in year

        Returns
        -------
        A dict of depreciation cost.

        '''

        i = self.interest
        operating_hours = self.operating_hours
        dry_mass = self._throughput
        depreciation_cost = {}
        for unit in self._units:
            for eqpt, cost in unit.purchase_costs.items():
                unit_eqpt = f'{unit.ID} - {eqpt}'
                n = unit.equipment_lifetime[eqpt]/(365*24)*operating_hours
                sal = unit.salvage[eqpt]
                frac = i*(1+i)**n/((1+i)**n-1)
                depreciation_cost[unit_eqpt] = \
                    (cost*(1-sal)*frac+cost*sal*i)/(operating_hours*dry_mass/1e3)
        return depreciation_cost

    def get_IHT_cost(self):
        r'''
        Calculate insurance, housing, and tax costs per metric tonne (MT) of dry feedstock.

        .. math:: cost[\frac{$}{yr}] = IHT_{frac} * \frac{P+S}{2}
        .. math:: cost[\frac{$}{MT}] = \frac{cost[\frac{$}{yr}]}{MT \ feedstock \ per \ yr}

        where:
            P: purchase cost

            S: salvage value = salvage fraction * purchase cost

            :math:`IHT_{frac}`: insurance, housing, and tax costs as a fraction of purchase cost

                :math:`IHT_{frac} = I_{frac} + H_{frac} + T_{frac}`

        Returns
        -------
        A dict of IHT cost.

        '''

        frac = self.insurance + self.housing + self.tax
        operating_hours = self.operating_hours
        dry_mass = self._throughput
        IHT_cost = {}
        for unit in self._units:
            for eqpt, cost in unit.purchase_costs.items():
                unit_eqpt = f'{unit.ID} - {eqpt}'
                IHT_cost[unit_eqpt] = \
                    cost*(1+unit.salvage[eqpt])/2*frac/(operating_hours*dry_mass/1e3)
        return IHT_cost

    def get_maintenance_cost(self):
        r'''
        Calculate maintenance cost per metric tonne (MT) of dry feedstock.

        .. math:: cost[\frac{$}{hr}] = \frac{P*maintenance_{frac}}{lifetime}
        .. math:: cost[\frac{$}{MT}] = \frac{cost[\frac{$}{hr}]}{dry \ feedstock \ per \ hour}


        where:
            P: purchase cost

            :math:`maintenance_{frac}`: maintenance cost as a fraction of purchase cost

            lifetime: equipment lifetime in hours

        Returns
        -------
        A dict of maintenance cost.

        '''

        dry_mass = self._throughput
        maintenance_cost = {}
        hr = self.operating_hours
        for unit in self._units:
            for eqpt, cost in unit.purchase_costs.items():
                unit_eqpt = f'{unit.ID} - {eqpt}'
                lifetime = unit.equipment_lifetime[eqpt]*hr
                maintenance_cost[unit_eqpt] = \
                    unit.maintenance[eqpt]*cost/lifetime/(dry_mass/1e3)
        return maintenance_cost

    def get_labor_cost(self):
        r'''
        Calculate labor cost per metric tonne (MT) of dry feedstock.

        .. math:: cost[\frac{$}{MT}] = \frac{cost[\frac{$}{hr}]}{dry \ feedstock \ per \ hour} * labor_adjustment


        Returns
        -------
        A dict of labor cost.

        '''

        adj = self.labor_adjustment
        labor_cost = {}
        for unit in self._units:
            for eqpt, cost in unit.purchase_costs.items():
                unit_eqpt = f'{unit.ID} - {eqpt}'
                if cost == 0:
                    labor_cost[unit_eqpt] = 0
                else:
                    labor_cost[unit_eqpt] = \
                        unit.labor[eqpt]/(unit.cost_items[eqpt].S/1e3)*adj
        return labor_cost

    def get_electricity_cost(self):
        '''
        Calculate electricity cost per metric tonne (MT) of dry feedstock.

        Returns
        -------
        A dict of electricity cost.

        '''

        e_dct = {}
        e_price = bst.PowerUtility.price
        for unit in self._units:
            for eqpt, cost_item in unit.cost_items.items():
                unit_eqpt = f'{unit.ID} - {eqpt}'
                if 'dryer' in eqpt.lower():
                    e_dct[unit_eqpt] = \
                        e_price * unit.power_utility.rate/(self._throughput/1e3)
                else:
                    e_dct[unit_eqpt] = e_price * cost_item.kW/(cost_item.S/1e3)
        return e_dct


    def get_equipment_cost(self):
        '''
        Get a result table for equipment costs, all in $ per metric
        tonne ([$/MT]), equivalent to [$/Mg].

        Returns
        -------
        A pandas.DataFrame table.

        '''
        cost_dct = {
            'Depreciation': self.get_depreciation_cost(),
            'IHT': self.get_IHT_cost(),
            'Maintenance': self.get_maintenance_cost(),
            'Labor': self.get_labor_cost(),
            'Electricity': self.get_electricity_cost()
            }
        df = pd.DataFrame.from_dict(cost_dct)
        df['Total'] = df.sum(axis='columns')
        all_series = df.sum(axis='index')
        all_series.name = 'All'
        df = df.append(all_series)
        return df

    def get_chemical_cost(self):
        '''
        Get a result table for chemical costs, all in $ per metric
        tonne ([$/MT]), equivalent to [$/Mg].

        Returns
        -------
        A pandas.DataFrame table.

        '''
        IDs = tuple(i.ID for i in self._feeds)
        mass_dct = dict.fromkeys(IDs)
        price_dct = dict.fromkeys(IDs)
        cost_dct = dict.fromkeys(IDs)
        dry_mass = self._throughput
        for feed in self._feeds:
            ID = feed.ID
            m = mass_dct[ID] = feed.F_mass/(dry_mass/1e3)
            p = price_dct[ID] = feed.price
            cost_dct[ID] = m * p
        chemical_dct = {
            'Mass [kg/MT]': mass_dct,
            'Price [$/kg]': price_dct,
            'Cost [$/MT]': cost_dct
            }
        df = pd.DataFrame.from_dict(chemical_dct)
        all_series = pd.DataFrame(columns=df.columns, index=('All',))
        all_series['Cost [$/MT]'] = df.sum(axis='index')['Cost [$/MT]']
        df = df.append(all_series)
        return df

    @property
    def depot_sys(self):
        '''
        [biosteam.System] Depot system.

        Note
        ----
        Changing the depot system will not automatically update feedstock.'''
        return self._depot_sys
    @depot_sys.setter
    def depot_sys(self, sys):
        self._depot_sys = sys
        self._units = sorted(sys.units, key=lambda u: u.ID)
        self._feeds = sorted((i for i in sys.feeds if i.price), key=lambda f: f.ID)

    @property
    def feedstock(self):
        '''[biosteam.Stream] Preporcessed feedstock stream from the depot system.'''
        return self._feedstock
    @feedstock.setter
    def feedstock(self, feed):
        self._feedstock = feed
        self._throughput = feed.F_mass - feed.imass['Water']

    @property
    def feedstock_unit_price(self):
        '''Feedstock price in $/metric tonne (MT), equivalent to 1 Mg.'''
        price0 = self.get_equipment_cost()['Total']['All']
        price1 = self.get_chemical_cost()['Cost [$/MT]']['All']
        return price0+price1



# %%

from biorefineries.ethanol_adipic._chemicals import chems
from biorefineries.ethanol_adipic._settings import _labor_2011to2016, price
bst.settings.set_thermo(chems)

flowsheet1 = create_default_depot(kind='CPP', with_AFEX=False)
flowsheet2 = create_default_depot(kind='HMPP', with_AFEX=False)
flowsheet3 = create_default_depot(kind='CPP', with_AFEX=True,
                            water_price=price['Makeup water'],
                            ammonia_price=price['NH3'],
                            natural_gas_price=price['Natural gas'])
flowsheet4 = create_default_depot(kind='HMPP', with_AFEX=True,
                            water_price=price['Makeup water'],
                            ammonia_price=price['NH3'],
                            natural_gas_price=price['Natural gas'])



flowsheet1.system.prep_sys.simulate()
flowsheet2.system.prep_sys.simulate()
flowsheet3.system.prep_sys.simulate()
flowsheet4.system.prep_sys.simulate()


cost1 = PreprocessingCost(depot_sys=flowsheet1.system.prep_sys,
                          labor_adjustment=_labor_2011to2016)
cost2 = PreprocessingCost(depot_sys=flowsheet2.system.prep_sys,
                          labor_adjustment=_labor_2011to2016)
cost3 = PreprocessingCost(depot_sys=flowsheet3.system.prep_sys,
                          labor_adjustment=_labor_2011to2016)
cost4 = PreprocessingCost(depot_sys=flowsheet4.system.prep_sys,
                          labor_adjustment=_labor_2011to2016)