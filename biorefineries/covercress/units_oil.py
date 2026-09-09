#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:24:13 2026

@author: princyk2
"""
"""
Reference
[1].https://acsess.onlinelibrary.wiley.com/doi/10.1002/csc2.21183#:~:text=The%20mean%20phenotypic%20responses%20averaged,for%20their%20effective%20practical%20utility.
[2]. Extraction of pennycress (Thlaspi arvense L.) seed oil by full pressing
[3]:Extraction of proteins from pennycress seeds and press cake
[4]Techno-economic analysis of hydroprocessed renewable jet fuel production
from pennycress oilseed
[5]https://www.oil-press-group.com/solvent-recovery-systems.html
[6].https://onecpm.com/success-story/the-complete-guide-to-canola-oil-processing-from-seed-to-oil
[7].Product and Process Design Principles: Synthesis, Analysis and Evaluation, 4th Edition
[8].https://onecpm.com/product/desolventizer-toaster---- 
[9]https://onecpm.com/product/vacuum-stripperfor stripper design
Warren D. Seider, Daniel R. Lewin, J. D. Seader, Soemantri Widagdo, Rafiqul Gani, Ka Ming Ng 
[10]https://www.alibaba.com/showroom/canola-seed-machine.html for euippemnt cost
[11]https://frenchoil.com/products/oilseed-equipment/
[12]https://www.aocs.org/resource/meal-desolventizing-toasting-drying-and-cooling/?
[13]https://www.feedingredientsasia.com/en/market-insights/applications-and-buyers/soy-protein-isolate-market-growth-2026?utm_source
[14].Leading Edge Technologies and Perspectives in Industrial Oilseed Extraction 
[15].https://junyuanpetroleumgroup.com/food-grade-hexane/oil-production-by-hexane-solvent-extraction/
[16].chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://pubs.rsc.org/en/content/articlepdf/1988/p2/p29880000523
[17].chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.nist.gov/system/files/documents/srd/jpcrd367.pdf
[18]. The Effect of Optimizing the Stripping and Drying Parameters During Industrial Extraction on the Physicochemical Properties of Soybean Oil
[19].Oil Press Group. (n.d.). Miscella evaporator. Retrieved May 26, 2026, from https://www.oil-press-group.com/miscella-evaporator.html
[20].https://onlinelibrary.wiley.com/doi/chapter-epub/10.1002/9781118535202.ch4
[21]https://www.klmtechgroup.com/PDF/Articles/articles/Steam%20Stripping%20Paper%20Version%20Final.pdf
[22].Measurements of activity coefficients at infinite dilution in vegetable oils and capric acid using the dilutor technique
[23]https://www.desmet.com/en/feed/solvent-extraction/distillation-and-solvent-recovery/sieve-tray-oil-stripper
[24]https://aocs.onlinelibrary.wiley.com/doi/full/10.1002/sfp2.1029
[25].Soybean oil extraction with ethanol from multiple-batch assays to reproduce a continuous, countercurrent, and multistage equipment
[26].Extraction of proteins from pennycress seeds and press cake
[27]"Cold-pressing, ethanol defatting, and saline extraction enhances properties of protein products from new pennycress varieties (covercress™)"
[28]Impact of extraction conditions and seed variety on the characteristics of pennycress (Thlaspi arvense) protein: a structure and function approach
[29].Refining Vegetable Oils: Chemical and Physical Refining
[30].Process simulation and techno-economic analysisof bio-jet fuel and green diesel production —Minimum selling prices
[31].Process modeling of hydrodeoxygenation to produce renewable jet fuel
and other hydrocarbon fuels
[32]Thermodynamic Equilibrium Analysis of Triolein Hydrodeoxygenation for
Green Diesel Production
[33].Process modeling of hydrodeoxygenation to produce renewable jet fuel
and other hydrocarbon fuels
[34].Technoeconomic analysis of biojet fuel production from camelina at
commercial scale: Case of Canadian Prairies
[35]. https://www.chemengonline.com/plant-cost-index-beta/    CEPCI index for 2026 & 2025
[36].Uncertainties in Early-Stage Capital Cost Estimation of Process Design – A Case Study on Biorefinery Design
[37]. https://www.researchgate.net/publication/289844027_Solid-Liquid_Extraction
[38].https://www.scirp.org/journal/paperinformation?paperid=89731
[39]https://chem-casts.com/tools/property-calculator/pure-component/14465-68-0-   # for critical properties of TAG and also for density and molar volue for TAG at diff high temp 
[40].chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.iata.org/en/iata-repository/publications/economic-reports/global-outlook-for-air-transport-june-2026/     # SAF price refer table 9
[41].chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https:/www3.epa.gov/ttnecas1/regdata/IPs/Vegetable%20Oil_IP.pdf
[42].Enzyme assisted protein extraction from rapeseed, soybean, and microalgae meals


"""
#reference--corn biorefinery (year 2007)
import biosteam as bst
import thermosteam as tmo
from  biosteam import Unit, HeatUtility
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year
from biosteam import tank_factory
from biosteam.exceptions import lb_warning
from biosteam import units
from biosteam.units.splitting import Splitter
from biosteam.units.drying import DrumDryer
from math import exp, log, ceil
from thermosteam import separations 
from thermosteam.reaction import Reaction, ParallelReaction, SeriesReaction
from biosteam.units.design_tools import PressureVessel
from biorefineries.covercress.utils import CEPCI
# from biorefineries import covercress
# from biorefineries.covercress import chemicals_cc
# from biorefineries.covercress.chemicals_cc import chemicals
"""
# %% Seed reception and cleaning units for covercress oilseed handling
"""
#since we use the refernce year 2007 the conveying belt has decorator with cost for year 2013 
#convert cost for belt from 2013 to 2007
CE_2013 = 567
CE_2007 = CEPCI_by_year[2007]
cost_2013_belt = 813
cost_2007_belt = cost_2013_belt * (CE_2007 / CE_2013) 
@cost('Flow rate', CE=CE_2007, cost=cost_2007_belt, ub=2000, n=0.38,
      BM=1.61, N='Number of conveyors')
class ConveyingBelt(Unit):
    length = 40 #: ft
    height = 20 #: ft
    _N_ins = _N_outs = 1
    _minimum_flow = 120
    _units = {'Flow rate': 'ft^3/hr'}
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
    
    def _design(self):
        feed = self.ins[0]
        self.design_results['Flow rate'] = F_vol = feed.F_vol*35.315 # ft3/hr
        if F_vol < self._minimum_flow:
            lb_warning(self, 'Flow rate', F_vol, 'kg/hr', 3, self._minimum_flow)
        F_mass = feed.F_mass * 0.0006124 # kg/hr to lb/s
        self.add_power_utility(
            0.7457 * (0.00058 * F_mass**0.82 * self.length + self.height*0.00182*F_mass) # kW
        )
#the final SI unit for power_utility isin KWh  
#ea
cost_2013_screen=1010
# TODO update the base cost as per the covercress 
cost_2007_screen = cost_2013_screen * (CE_2007 / CE_2013)
@cost('Area', ub=200, CE=CE_2007, cost= cost_2007_screen , n=0.91, BM=1.73, N='Number of screens')
class VibratingScreen(Splitter):
    # Assume 3-deck vibrating screen
    # Mean seed dimensions for pennycress phenotypes (mm): length 1.34, width 1.79 refer [1]
    seed_length_mm = 1.34 #TODO changed the seed length and width, as per the covercress seed dimension 
    seed_width_mm = 1.79 
    # Flow rate per area of screen per aperture (kg/(ft^2-hr-mm))
    capacity = 0.2 * 907.18474
    #0.2 is (0.2 is the vendor‑supplied capacity constant)capacity in US SI units (lb) from a correlations from vendor
    #1 short ton = 2000 lb
    #1 lb = 0.45359237 kg
    #1 short ton = 2000 × 0.45359237 ≈ 907.18474 kg
    # Mesh opening (mm): set slightly larger than seed so seeds pass, debris retained
    # e.g. ~1.4× larger dimension so opening ≈ 2.5 mm #TODO the meshopening
    mesh_opening = 2.5 
    _units = {'Area': 'ft^2'}
    def _design(self):
        feed = self.ins[0]
        self.design_results['Area'] = feed.F_mass / (self.capacity * self.mesh_opening)
        #the cost decorator use the design result variable ie Area which then check the ub if its >200 then (its the maximum area for one unit) then it choose for multiple units by calc number of screen thus the total purhcase cost reflects the number of parralel units.
cost_2013_mag=1010  
cost_2007_mag = cost_2013_mag * (CE_2007 / CE_2013)
@cost('Flow rate', units='kg/hr', CE=CE_2007 , cost= cost_2007_mag , S=333333, BM=4.16, n=0.6)
class MagneticSeparator(Unit): 
    _N_outs = 2
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
        
@cost('Flow rate', units='kg/hr', CE=CE_2007, cost=60300., S=45350., n=0.6, ub=7.2e5, BM=4)
class SeedCleaningSystem(Splitter):
    """
    Cleaning system (aspiration / fines removal). One feed, two outs: [0] clean seed, [1] fines/dust.
    Splits (e.g. isplit) must be set on the flowsheet so the right components go to each outlet.
    """
    _N_ins = 1
    _N_outs = 2
    _units = {'Flow rate': 'kg/hr'}
    def _design(self):
        feed = self.ins[0]
        self.design_results['Flow rate'] = feed.F_mass
SeedStorage = tank_factory(
    'SeedStorage',
    CE=CE_2007,
    cost=979300.,      # reference purchase cost at S
    S=185400.,        # reference size (effective volume or flow×tau)
    tau=259.2,        # residence time (h)
    n=1.0,
    V_wf=0.9,
    V_max=3e5,        # max single-tank volume (m³)
    V_units='m3',
    BM=4,
)  #the parameters imported from the corn ->corn storage

#%% Seed preparation
"""
#%% Seed preparation
"""


#TODO update  refer [2] for cost ,later the cost-replace with https://www.alibaba.com/showroom/oil-seed-cooker.html change the cost and corresponding area S
@cost('Peripheral drum area', ID='Conditioning_Drum', CE=CEPCI_by_year[2007], ub=7854.0, BM=2.06,
      S=1235.35, units='m2', n=0.6, cost=0.52 * 2268000., kW=938.866)
class SeedConditioning(DrumDryer): #conditioning unit
      cost_items = {}   # clear parent DrumDryer purchase-cost entries
      pass

       
        

# TODO: Replace `cost` / `S` / `kJ_per_kg` with vendor data or literature for your roller flaker. asssumed cost of roll presse from the seider 
@cost('Flow rate', units='lb/hr', CE=567, lb=150, ub=12000, BM=1.39, 
     f=lambda S: exp(10.9807 - 0.4467*log(S) + 0.06136*log(S)**2))
class SeedFlaking(Unit): #defiened flaking calss using roller press cost correlation adn as like screw press inherited form the solids separtor
    _N_ins=1
    _N_outs=1
    
    kWh_per_bmt = 37.2 # From Perry's Handbook, 18-126 also equal to 
    kJ_per_kg=kWh_per_bmt *3.6 # converted KWh to KJ by 37.2*3600/1000
    _units = {'Flow rate': 'lb/hr'}
    def _design(self):
        feed, = self.ins
        self.design_results['Flow rate'] = feed.F_mass
        if feed.F_mass<=0:
           return
        power_kW = self.kJ_per_kg * feed.F_mass / 3600.0
        self.add_power_utility(power_kW)
       
    
    def _run(self):
        feed, =self.ins
        flake_out,=self.outs
        if feed.F_mass<=0:
           flake_out.copy_like(feed)
           return
        flake_out.copy_like(feed)
        # kW = (kJ/kg) × (kg/hr) / (3600 s/hr)
    def _cost(self):
        self._decorated_cost()
        biomass = self.ins[0]
        bmt = biomass.F_mass * 0.001 #kg-metrictonnes
        self.add_power_utility(bmt * self.kWh_per_bmt)

        
#TODO update  refer [2] for cost ,later the cost-replace with https://www.alibaba.com/showroom/oil-seed-cooker.html change the cost and corresponding area S
@cost('Peripheral drum area',ID='Cooking_Drum', CE=CEPCI_by_year[2007], ub=7854.0, BM=2.06,
      S=1235.35, units='m2', n=0.6, cost=0.52 * 2268000., kW=938.866)      
class Cooking(DrumDryer): #cooking unit refer[2]
      cost_items = {}   # clear parent DrumDryer purchase-cost entries
      pass
# @cost('Flow rate', 'Pre-press',units='kg/hr',cost=42_000,S=8333, ub=8333,n=0.6,kW=30,BM=1.65,CE=858,) 
    #TODO update CEPCI for 2026
    
 #https://zzjtian.en.made-in-china.com/product/dXkmQKAVHepr/China-Palm-Oil-Residue-Spent-Grain-Dewatering-Screw-Filter-Press-for-Sale.html?pv_id=1jrbd2pnk219&faw_id=1jrbdacmp6ce&bv_id=1jrbecckm8d4&pbv_id=1jrbd2o3r2a5   
             # reference: max is 200 t/day# 
             # max per machine: 200 t/day cost decorator has line ceil which calcaute how many number of screw press is required because the quote says max is 200t/day but this model has cacapity >200t/day so we need more than 1 screwpress
    #BM=1.65 for crushing 
    
    
    

# class CovercressScrewPress(bst.ScrewPress):
#       # cost_items = {}
#       # pass 
#       def _design(self):
#           feed, = self.ins
#           self.design_results['Flow rate'] = feed.F_mass  # kg/hr  compared to ub

#       def _cost(self):
#         self._decorated_cost()   #  parallel logic runs here solidsepartor has not method to calcaute the cost 
# TODO: Replace cost / S / kW with vendor quotes or hexane-specific references.updated from seider check its correct
#vapour pressure of hexane is 0.26psig
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=29900., S=10000.0, CE=CE_2007, n=0.513, BM=1) # refer[seider]page no 485 for conerrof for carbon steel  where V =10,0000 n= 0.513, Bm=1 for carbonsteel refer pgno 480
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=7493, S=1981, CE=CEPCI[2010], n=0.8, BM=2.3) #TODO update: refered from sulphuric storage tank succinic_units
class HexaneStorageTank(Unit):
    """
    Liquid hexane storage with tank + transfer pump capital (flow-rate basis).
    Inlet and outlet streams are the same composition; use for makeup or recycle lines.
    """

    _N_ins = _N_outs = 1
    _units = {'Flow rate': 'kg/hr'}

    def _run(self):
        self.outs[0].copy_like(self.ins[0])

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=29900., S=10000.0, CE=CE_2007, n=0.513, BM=1) # refer[seider]page no 485 for conerrof for carbon steel  where V =10,0000 n= 0.513, Bm=1 for carbonsteel refer pgno 480
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=7493, S=1981, CE=CEPCI[2010], n=0.8, BM=2.3) #TODO update: refered from sulphuric storage tank succinic_units
class EthanolStorageTank(Unit):
    """
    Liquid hexane storage with tank + transfer pump capital (flow-rate basis).
    Inlet and outlet streams are the same composition; use for makeup or recycle lines.
    """

    _N_ins = _N_outs = 1
    _units = {'Flow rate': 'kg/hr'}

    def _run(self):
        self.outs[0].copy_like(self.ins[0])

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass
    
 # Extraction of oil from cake refer[3]  and refer[4]
solids_retentate={'Glucose','Protein','Ash','Lignin','Sucrose','Cellulose','Hemicellulose'}
lipids = (
     'Lipid', 'Oil',
     'PL', 'Phospholipid',
     'FFA', 'OleicAcid',
     'MAG', 'MonoOlein', 'DAG', 'DiOlein', 'TAG', 'TriOlein',
     'LinoleicAcid', 'LinolenicAcid','Stericacid','OOO', 'LLL', 'LnLnLn','SSS',
 )

# Rotocel or Carousel extractor the exctractor used in the plant refer [37]fig 6.6
#TODO cost decorator check for any auxallaries are utilized by the industrial grade countercurrent multistage extractor
#TODO later the process should also compre with multistagemixersettler or solid liquid liquid equilibrium centrifuge (collect partition coefficent data )
@cost('Flow rate', units='kg/hr', CE=858, cost=180_000, S=3333,ub=3333, n=0.6, BM=4.16) # 2026 cost from above refer BM =4.16 assumed for vertical pressure vessels from refer waren sider table 16.11
 #https://chinaoilmill.en.made-in-china.com/product/cvonaFRYJJrb/China-Rotocel-Castor-Seed-Extraction-Equipment-Soya-Bean-Cake-Edible-Oil-Solvent-Extractor.html?pv_id=1jrdsfub1e55&faw_id=1jrdskp0fffb&bv_id=1jrdskp0gb24&pbv_id=1jrdsfst055ev
 #  why n=0.6 refer[36] 0.60 is the standard default, the typical working range 
class CovercressExtractor(bst.units.LLEUnit):
    _N_ins = 1
    _N_outs = 2  # [0]=miscella (light), [1]=cake (heavy)

    def _init(self, top_chemical='Hexane', efficiency=1.0, **kwargs):
        forced_split_IDs = (
            'Cellulose', 'Hemicellulose', 'Lignin', 'Ash', 'Protein',
            'Glucose', 'Sucrose',
            'Hexane', 'OOO', 'LLL', 'LnLnLn','SSS',
            'Water'
        )
        forced_split = (
            0.0, 0.0, 0.0, 0.0, 0.0,  # solids to cake
            0.0, 0.0,                      # sugars to cake
            0.2, 0.999, 0.999, 0.999,0.999,        # solvent/oils to miscella
            0.049,                          # water carryover to miscella
        )
        super()._init(
            top_chemical=top_chemical,
            efficiency=efficiency,
            forced_split_IDs=forced_split_IDs,
            forced_split=forced_split,
            **kwargs,
        )
    def _design(self):
        feed, = self.ins
        self.design_results['Flow rate'] = feed.F_mass   # basis for @cost correlation
# Because true thermodynamic equilibrium represents a theoretical limit that requires an infinite amount of time or an infinite number of ideal stages to achieve
# . Real industrial decanters and extraction columns are imperfect. Due to fluid mechanics, mixing limitations, and emulsion formations, droplets of one liquid will often remain physically entrained in the other. This code applies a stage efficiency correction
# , forcing the software to mimic the imperfect, non-equilibrium mechanical separation of real equipment while still strictly obeying the Conservation of Mass.

#solvent recovery refer[5]
#TODO update cost for DT
@cost('Peripheral drum area',ID='Desolventizer_Drum', CE=CEPCI_by_year[2007], ub=7854.0, BM=2.06,
      S=1235.35, units='m2', n=0.6, cost=0.52 * 2268000., kW=938.866)      
class Desolventizer(DrumDryer): #cooking unit refer[2]
      cost_items = {}   # clear parent DrumDryer purchase-cost entries
      pass
#TODO add cost of drying and coo,ling unit
@cost('Peripheral drum area',ID='dryer_Drum', CE=CEPCI_by_year[2007], ub=7854.0, BM=2.06,
       S=1235.35, units='m2', n=0.6, cost=0.52 * 2268000., kW=938.866)      
class Drying(DrumDryer): #cooking unit refer[2]
       cost_items = {}   # clear parent DrumDryer purchase-cost entries
       pass 
class MealCooler(bst.HXutility):
    line = 'Meal cooler'

    # def _cost(self):
    #     """Skip HX purchase correlation DT/DC cost is on F102 Drying.""" 
    #     return
class DTDC(bst.RotaryVacuumFilter):
           pass
##
#protein recovery units
###


# Default kW based on an industrial corn dry-grind hammer mill
@cost('Flow rate', units='ton/hr', cost=4310, lb=2, ub=200,
      CE=567, n=0.78, kW=6.17, BM=2.3) 
class HammerMill(Unit):  pass

@cost(basis='Flow rate', ID='Ultrafilter', units='m3/hr',
      cost=2048000*0.0297, S=1303*0.2271246, CE=567, n=0.6, BM=2.5) #refer units, ethnaol_adipic muconic membrane
class ProteinUltrafiltration(bst.Unit): # refer[28] during ultrafiltraion feed is in liquid we cant use rortary vaccum filter becuase the feed for RVF is solid
    _N_ins = 1
    _N_outs = 2   # [0] permeate, [1] retentate (protein concentrate)

    def _init(self, protein_retention=0.99, salt_retention=0.02):
        self.protein_retention = protein_retention
        self.salt_retention = salt_retention

    def _run(self):
        feed, = self.ins
        permeate, retentate = self.outs
        permeate.empty()
        retentate.empty()
       

        # protein -> retentate 
        P = feed.imass['Protein']
        retentate.imass['Protein'] = self.protein_retention * P
        permeate.imass['Protein']  = (1 - self.protein_retention) * P

        # salts mostly -> permeate
        for ID in ('NaCl', 'NaOH'):
            if ID in feed.chemicals:
                retentate.imass[ID] = self.salt_retention * feed.imass[ID]
                permeate.imass[ID]  = (1 - self.salt_retention) * feed.imass[ID]

        # small solubles -> mostly permeate
        for ID in ('Glucose', 'Sucrose', 'Hexane'):
            if ID in feed.chemicals:
                retentate.imass[ID] = 0.05 * feed.imass[ID]
                permeate.imass[ID]  = 0.95 * feed.imass[ID]

        # water -> mostly permeate
        retentate.imass['Water'] = 0.05 * feed.imass['Water']
        permeate.imass['Water']  = 0.95 * feed.imass['Water']

        # anything else
        accounted = {'Protein','NaCl','NaOH','Glucose','Sucrose','Hexane','Water'}
        for ID in feed.chemicals.IDs:
            if ID not in accounted:
                permeate.imass[ID] += feed.imass[ID]

        retentate.T = permeate.T = feed.T

    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_vol
from biosteam.units.drying import DrumDryer  # only for sizing helpers if you want
from biosteam.units.design_tools import cylinder_area, cylinder_diameter_from_volume
from biosteam import HeatUtility

_chilled = HeatUtility.get_cooling_agent('chilled_water')  # same as process_settings

@cost('Peripheral drum area', ID='Freeze_dryer_drum', CE=CEPCI_by_year[2007],
      ub=7854.0, BM=2.06, S=1235.35, units='m2', n=0.6,
      cost=0.52 * 2268000., kW=200.)   # calibrate kW to vacuum + refrigeration
class Lyophilizer(Unit):   # used drumdryr= cost and design and run to calcu uses separions.
    _N_ins = 1
    _N_outs = 2
    _units = {'Peripheral drum area': 'm2'}

    def _init(self, moisture_content=0.001, T=233.15, H=20., length_to_diameter=20):
        # H and Leng to dia same as drum dryer
        self.moisture_content = moisture_content
        self.T = T
        self.H = H
        self.length_to_diameter = length_to_diameter
        self.h_fg_sublimation = 2830.  # kJ/kg water removed (lit. value)

    def _run(self):
        feed, = self.ins
        powder, vapor = self.outs

        powder.copy_like(feed)
        vapor.empty() # make sure no water vapour present before run the separtions to keep mass balance
        vapor.phase = 'g'

        # Move excess water from powder -> vapor until target moisture
        separations.adjust_moisture_content(powder, vapor, self.moisture_content) # used separation module to  enforces a target water fraction by moving water between two outlet streams.

        powder.T = vapor.T = self.T   # low-temperature operation

    def _design(self):
        water_removed = self.outs[1].imass['Water']  # kg/hr
        self.design_results['Evaporation'] = water_removed
        volume = water_removed / self.H if self.H > 0 else 0
        D = cylinder_diameter_from_volume(volume, self.length_to_diameter)
        L = D * self.length_to_diameter
        self.design_results['Peripheral drum area'] = cylinder_area(D, L)

        # Cooling duty: sublimation + sensible (use chilled water, not steam)
        if water_removed > 0:
            duty = water_removed * self.h_fg_sublimation  # kJ/hr
            self.add_heat_utility(duty, T_in=self.T, T_out=self.T, agent=_chilled)
            # optional: vacuum pump electricity
            self.add_power_utility(0.5 * water_removed / 3600)  # tune kWh/kg        
# @cost(basis='Flow rate', ID='Lyophilizer', units='m3/hr',
#       CE=CE_2007, cost=500000., S=50., n=0.6, BM=2.5) #TODO update the cost 
# class Lyophilizer(Unit):
#     """
#     Freeze-drying proxy for UF retentate.
#     Removes water to very low moisture at low T.
#     Capital: Lyophilizer purchase cost.
#     Operating: electricity from kWh per kg water removed.
#     """
#     line = 'Lyophilizer'
#     _N_ins = 1
#     _N_outs = 2   # [0] protein powder, [1] water vapor

#     # Electricity intensity [kWh/kg water removed] — #TODO update from literature/vendor
#     kWh_per_kg_water = 25.0

#     _units = {'Flow rate': 'm3/hr'}

#     def _init(self, moisture_content=0.03, T=233.15):
#         self.moisture_content = moisture_content   # very dry powder (3 wt% water)
#         self.T = T                                 # low shelf temperature [K]

#     def _run(self):
#         feed, = self.ins
#         powder, vapor = self.outs

#         powder.copy_like(feed)
#         vapor.empty() # make sure no water vapour present before run the separtions to keep mass balance
#         vapor.phase = 'g'

#         # Move excess water from powder -> vapor until target moisture
#         separations.adjust_moisture_content(powder, vapor, self.moisture_content)

#         powder.T = vapor.T = self.T   # low-temperature operation

#     def _design(self):
#         feed, = self.ins
#         self.design_results['Flow rate'] = feed.F_vol   # m3/hr — capital cost basis

#         water_removed = self.outs[1].imass['Water']     # kg/hr
#         if water_removed > 0:
#             self.add_power_utility(self.kWh_per_kg_water * water_removed / 3600)  # kWh-kW
#     def _design(self):
#         length_to_diameter = self.length_to_diameter
#         design_results = self.design_results
#         design_results['Volume'] = volume = design_results['Evaporation'] / self.H 
#         design_results['Diameter'] = diameter = cylinder_diameter_from_volume(volume, length_to_diameter)
#         design_results['Length'] = length = diameter * length_to_diameter
#         design_results['Peripheral drum area'] = cylinder_area(diameter, length)
#         if self.utility_agent == 'Steam':
#             self.add_heat_utility(self.H_out - self.H_in, self.T)
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=196000, S=1171/17.031*35.046, CE=CEPCI[2010], n=0.7, BM=2)
class HydrogenGasStorageTank(bst.StorageTank):
    """H2 storage tank (BDO costing on BioSTEAM StorageTank)."""
    pass
TAG_IDs = ('OOO', 'LLL', 'LnLnLn', 'SSS')           
# @cost('Flow rate', 'Reactor',CE=567, cost=50_000, S=10_000, n=0.7, kW=1.5, BM=4.16,) #TODO update later  assumtion 0. pressure vessel refer diag fig 16.13 and pg 465 (autoclave cost can be used for F_vol till 4000 gal our capcaity is > 25106 so we decided to use cost of pressure vessel1. autoclave cost from the plant economic book Sieder 481 and ~ 1222 gal/hr and Vertical pressurevessels 4.16 table 16.11 from process economics 2. cost from refer [camelina] and scalin exponent for pressure vessel from humbird table 26 , pg 60

from math import pi, ceil
from biosteam.units.design_tools import PressureVessel
from biosteam.exceptions import DesignError

# refer HP the design of pressure vessel both reactors are at high pressure so choose design for pressure vessel  
class Hydrodeoxygenation(bst.Unit, PressureVessel):
    _N_ins = 3
    _N_outs = 1
    _units = {**PressureVessel._units, 'Total volume': 'm3', 'Number of reactors': ''}
    _V_max = pi/4*(20**2)*40/35.3147   # V_masx is the maximum volume of a reactor from eq= V=pi/4 D^2L where D=20 ft max dia and L=40 ft max lenth refer HP
    def __init__(self, ID='', ins=None, outs=(), T=273+310, P=3e6, tau=2.0,
             catalyst_wt_frac=0.05, excess_H2=0, rxns=None,
             V_wf=0.8, length_to_diameter=2, kW_per_m3=0.985,
             wall_thickness_factor=1,
             vessel_material='Stainless steel 316', vessel_type='Vertical'):
        Unit.__init__(self, ID, ins, outs)
        self.T, self.P, self.tau = T, P, tau
        self.catalyst_wt_frac = catalyst_wt_frac
        self.excess_H2 = excess_H2
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor   # <-- ADD THIS
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
            
# class Hydrodeoxygenation(bst.Unit):
#     _bounds={'Volume': (0.1,50)}
#     _units={'Volume':'m^3'}
#     _N_ins=3
#     _N_outs=1 # refer [30] PFD where outlet is 1, the outs then undergoes separions
#     def __init__(self,ID='', ins=None, outs=(), T=273+310, P=3e6, tau=2.0, catalyst_wt_frac=0.05,excess_H2=0, rxns=None,): #refer[4]
#         Unit.__init__(self,ID,ins,outs)
#         self.T=T
#         self.P=P
#         self.tau=tau
#         self.catalyst_wt_frac=catalyst_wt_frac
#         self.excess_H2=excess_H2
        if rxns is None:
            self.saturation_rxns = ParallelReaction([
                Reaction('OOO + 3 Hydrogen -> SSS', reactant='OOO', X=1.0),
                Reaction('LLL + 6 Hydrogen -> SSS', reactant='LLL', X=1.0),
                Reaction('LnLnLn + 9 Hydrogen -> SSS', reactant='LnLnLn', X=1.0),
            ])
            self.cracking_rxn = Reaction(
                'SSS + 3 Hydrogen -> 3 Stericacid + Propane', reactant='SSS', X=1.0
            )
            self.ffa_rxns = ParallelReaction([
                Reaction('Stericacid + 3 Hydrogen -> Octadecane + 2 Water', reactant='Stericacid', X=0.05),
                Reaction('Stericacid + Hydrogen -> Heptadecane + Water + CO', reactant='Stericacid', X=0.05),
                Reaction('Stericacid -> Heptadecane + CO2', reactant='Stericacid', X=0.90),
            ]) #parrallel because when several reactions shared the same reactant selectovity is different refer[4] for conversion rates
        else:
            self.rxns=rxns
    def _required_H2_mol(self,feed):
        ooo,lll,ln,sss=(feed.imol[i] for i in TAG_IDs)
        h2_for_saturation=3*ooo+ 6*lll +9*ln  #refer[ppt] #TODO check how moles calcuted in the reaction
        sss_total= ooo + lll+ln+sss
        h2_for_cracking =3*sss_total
        stericacid=3*sss_total
        h2_for_decarboxylation= stericacid*(3*0.05+0.05)
        return h2_for_saturation+h2_for_cracking+h2_for_decarboxylation
    def _required_catalyst(self,feed):
        oil_mass=sum(feed.imass[i] for i in TAG_IDs)
        return self.catalyst_wt_frac*oil_mass
    def _run(self):
        feed,h2,catalyst=self.ins
        product,=self.outs
        h2.phase='g'
        h2.imol['Hydrogen']=(1+self.excess_H2)*self._required_H2_mol(feed)
        # catalyst.imass['Pt/Al2O3'] = self._required_catalyst(feed)
        product.mix_from([feed,h2,catalyst],energy_balance=False) #. Real reactors are temperature-controlled, not adiabatic mixers we are not letting the mixture to calcute the etemop from the mixture enthalpy and also the reactor wont operate in the operting temp it will run at some intermediate  T 
        self.saturation_rxns(product)
        self.cracking_rxn(product)# product formed based on the reaction
        self.ffa_rxns(product)
        product.T=self.T
        product.P=self.P
    # def _design(self):
    #     effluent=self.outs[0]
    #     #assumed working volume as 80 %
    #     self.design_results['Flow rate'] = effluent.F_mass 
    #     # self.design_results['Volume']= self.tau*effluent.F_vol/0.8 
    #     self.add_heat_utility(self.Hnet,effluent.T)
    def _design(self):
        Design = self.design_results
        ins_F_vol = self.ins[0].F_vol   # preheated liquid oil only (HX301-1) or if we write F_vol_in with calactlsyt biosteam crash as it fails to find liquid molar enthalpy of catalsyt
        V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor
        
        N = ceil(V_total/self._V_max)
        if N == 0:
            V_reactor = 0
            D = 0
            L = 0
        else:
            V_reactor = V_total / N
            D = (4*V_reactor/pi/length_to_diameter)**(1/3)
            D *= 3.28084 # convert from m to ft
            L = D * length_to_diameter

        Design['Residence time'] = self.tau
        Design['Total volume'] = V_total
        Design['Single reactor volume'] = V_reactor
        Design['Number of reactors'] = N
        P, D, L = float(P), float(D), float(L)
        Design.update(self._vessel_design(P, D, L))
        if wall_thickness_factor == 1: pass
        elif wall_thickness_factor < 1:
            raise DesignError('wall_thickness_factor must be larger than 1')
        else:
              Design['Wall thickness'] *= wall_thickness_factor
              # Weight is proportional to wall thickness in PressureVessel design
              Design['Weight'] = round(Design['Weight']*wall_thickness_factor,2)
            
    def _cost(self):
        Design = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        
        if Design['Total volume'] == 0:
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] = 0
        
        else:
            baseline_purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] *= Design['Number of reactors']
            
            self.power_utility(self.kW_per_m3 * Design['Total volume'])

paraffin_IDs = ('Octadecane', 'Heptadecane')

# @cost('Volume', 'Reactor',CE=567, cost=2250, n=0.58, kW=1.5, BM=4.16,)#cost from refer [camelina] and scalin exponent for pressure vessel from humbird table 26 , pg 60 for 2 price TCI was very high so we used autoclave cost
class Hydrocracking(bst.Unit, PressureVessel):
    """
    Isomerisation-reactor hydrocracking (ref [4]).
    Eq. (9)  Cn + H2 -> C(n-3) + C3H8
    Eq. (10) Cn + H2 -> C(n-8) + C8H18
    Eq. (11) Cn + H2 -> C(n-7) + C7H16
    Eq. (12) C18 + H2 -> 2 C9H20
   
    """
    _N_ins = 3
    _N_outs = 1
    _units = {**PressureVessel._units, 'Total volume': 'm3', 'Number of reactors': ''}
    _V_max = pi/4*(20**2)*40/35.3147   # V_masx is the maximum volume of a reactor from eq= V=pi/4 D^2L where D=20 ft max dia and L=40 ft max lenth refer HP
    def __init__(self, ID='', ins=None, outs=(), T=273+310, P=3e6, tau=2.0,
             catalyst_wt_frac=0.05, excess_H2=0, rxns=None,
             V_wf=0.8, length_to_diameter=2, kW_per_m3=0.985,
             wall_thickness_factor=1,
             vessel_material='Stainless steel 316', vessel_type='Vertical'):
        Unit.__init__(self, ID, ins, outs)
        self.T, self.P, self.tau = T, P, tau
        self.catalyst_wt_frac = catalyst_wt_frac
        self.excess_H2 = excess_H2
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor   # <-- ADD THIS
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
    # _bounds = {'Volume': (0.1, 50)}
    # _units = {'Volume': 'm^3'}
    # _N_ins = 3   # feed, H2, catalyst makeup
    # _N_outs = 1

    # def __init__(self, ID='', ins=None, outs=(),
    #              T=273.15 + 250, P=3e6, tau=1.5,
    #              catalyst_wt_frac=0.02, excess_H2=0.3):
    #     Unit.__init__(self, ID, ins, outs)
    #     self.T = T
    #     self.P = P
    #     self.tau = tau
    #     self.catalyst_wt_frac = catalyst_wt_frac
    #     self.excess_H2 = excess_H2

        # Literature pathway conversions (molar % as fraction)
        # C17: eq9=50.4, eq10=17.3, eq11=2.0
        # C18: eq9=4.9,  eq10=41.9, eq11=13.9, eq12=23.3
        self.hck_rxns = ParallelReaction([
            Reaction('Octadecane + Hydrogen -> Pentadecane + Propane',
                     reactant='Octadecane', X=0.049),
            Reaction('Heptadecane + Hydrogen -> Tetradecane + Propane',
                     reactant='Heptadecane', X=0.504),

            Reaction('Octadecane + Hydrogen -> Decane + Octane',
                     reactant='Octadecane', X=0.419),
            Reaction('Heptadecane + Hydrogen -> Nonane + Octane',
                     reactant='Heptadecane', X=0.173),

            Reaction('Octadecane + Hydrogen -> Undecane + Heptane',
                     reactant='Octadecane', X=0.139),
            Reaction('Heptadecane + Hydrogen -> Decane + Heptane',
                     reactant='Heptadecane', X=0.020),

            Reaction('Octadecane + Hydrogen -> 2 Nonane',
                     reactant='Octadecane', X=0.233),
        ])
        self.isom_rxns = ParallelReaction([
    # n = 8–16 ([30] table 4) iso C7=0 same refernce
    Reaction('Octane -> IsoOctane',             reactant='Octane',       X=0.90),
    Reaction('Nonane -> IsoNonane',             reactant='Nonane',       X=0.88),
    Reaction('Decane -> IsoDecane',             reactant='Decane',       X=0.90),
    Reaction('Undecane -> IsoUndecane',         reactant='Undecane',     X=0.90),
    # Reaction('Dodecane -> IsoDodecane',         reactant='Dodecane',     X=0.068),
    # Reaction('Tridecane -> IsoTridecane',       reactant='Tridecane',    X=0.076),
    Reaction('Tetradecane -> IsoTetradecane',   reactant='Tetradecane',  X=0.539),
    Reaction('Pentadecane -> IsoPentadecane',   reactant='Pentadecane',  X=0.272),
   # Reaction('Hexadecane -> IsoHexadecane',     reactant='Hexadecane',   X=0.100),
    # n = 17–18 ( table 4 [30])
    Reaction('Heptadecane -> IsoHeptadecane',     reactant='Heptadecane',  X=0.126),
    Reaction('Octadecane -> IsoOctadecane',     reactant='Octadecane',   X=0.189),
])

    def _X_sum(self, reactant):
        return sum(r.X for r in self.hck_rxns if r.reactant == reactant)

    def _required_H2_mol(self, feed):
        n18 = feed.imol['Octadecane']
        n17 = feed.imol['Heptadecane']
        return n18 * self._X_sum('Octadecane') + n17 * self._X_sum('Heptadecane')

    def _required_catalyst(self, feed):
        paraffin_mass = sum(feed.imass[i] for i in paraffin_IDs)
        return self.catalyst_wt_frac * paraffin_mass

    def _run(self):
        feed, h2,catalyst = self.ins
        product, = self.outs
        
        h2.phase = 'g'
        h2.imol['Hydrogen'] =(1 + self.excess_H2) * self._required_H2_mol(feed) #
        catalyst.phase = 's'
        # catalyst.imass['Pt/Al2O3'] = self._required_catalyst(feed)
        product.mix_from([feed,h2,catalyst], energy_balance=False)
        # why energy_balane=false#. Real reactors are temperature-controlled, not adiabatic mixers we are not letting the mixture to calcute the etemop from the mixture enthalpy and also the reactor wont operate in the operting temp it will run at some intermediate  T 
        self.hck_rxns(product.mol)
        self.isom_rxns(product.mol) 
        product.T = self.T
        product.P = self.P
    def _design(self):
        Design = self.design_results
        ins_F_vol = self.ins[0].F_vol   # preheated liquid oil only (HX301-1) or if we write F_vol_in with calactlsyt biosteam crash as it fails to find liquid molar enthalpy of catalsyt
        V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor
        
        N = ceil(V_total/self._V_max)
        if N == 0:
            V_reactor = 0
            D = 0
            L = 0
        else:
            V_reactor = V_total / N
            D = (4*V_reactor/pi/length_to_diameter)**(1/3)
            D *= 3.28084 # convert from m to ft
            L = D * length_to_diameter

        Design['Residence time'] = self.tau
        Design['Total volume'] = V_total
        Design['Single reactor volume'] = V_reactor
        Design['Number of reactors'] = N
        P, D, L = float(P), float(D), float(L)
        Design.update(self._vessel_design(P, D, L))
        if wall_thickness_factor == 1: pass
        elif wall_thickness_factor < 1:
            raise DesignError('wall_thickness_factor must be larger than 1')
        else:
              Design['Wall thickness'] *= wall_thickness_factor
              # Weight is proportional to wall thickness in PressureVessel design
              Design['Weight'] = round(Design['Weight']*wall_thickness_factor,2)
                
    def _cost(self):
        Design = self.design_results
        baseline_purchase_costs = self.baseline_purchase_costs
        
        if Design['Total volume'] == 0:
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] = 0
        
        else:
            baseline_purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] *= Design['Number of reactors']
            
            self.power_utility(self.kW_per_m3 * Design['Total volume'])
    # def _design(self):
    #     effluent=self.outs[0]
    #     # don’t let gas/supercritical components drive liquid V in sizing — filter them out at design time.Only _design skips it for F_vol (avoids COSTALD at T > Tc).
    #     # thus avoiding the R302 crash when it fails to evalutte the liquid molar volume  for hexane at supercritical stage and for simialr gases so design ignore gases 
    #     # effluent = self.outs[0]
       
    #     # gases = ('Hexane', 'Hydrogen', 'Pt/Al2O3', 'CO2', 'CO', ) #'Water', 'Propane')
    #     # saved = {ID: effluent.imol[ID] for ID in gases if ID in effluent.chemicals.IDs}
    #     # try:
    #     #     for ID in saved:
    #     #         effluent.imol[ID] = 0
    #     #     F_vol = effluent.F_vol if effluent.F_mol else 0.
    #     # finally:
    #     #     for ID, n in saved.items():
    #     #         effluent.imol[ID] = n
    #     self.design_results['Volume'] = self.tau * effluent.F_vol / 0.8
    #     self.add_heat_utility(self.Hnet, effluent.T)
    #     effluent = self.outs[0]
    #     self.design_results['Volume'] = self.tau * effluent.F_vol / 0.8
    #     self.add_heat_utility(self.Hnet, effluent.T)    
    # # def _cost(self):
    # #     D = self.design_results
    # #     self.purchase_costs.update(
    # #         self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
    # #     )
            
            
            
@cost('Propane flow rate', S=1.0, CE=596.2, cost=5.0e6,   n=0.6,  kW=0.5 * 3600,)  # design basis = recovered C3, NOT H2                 # kg/s reference — set from your reference [34]
              # USD at CE — replace with cited value
    
  
            # distilaltion unit cost 

class PropaneDistillation(bst.Splitter):
    """
    Propane (LPG) recovery from H2/CO2/CO off-gas.
    Surrogate: compression, knock-out, and refrigeration to condense C3.
    Split fractions from refer [34]; capital scaled on propane product flow.
    """
    _units = {'Propane flow rate': 'kg/s'}
    _line = 'Propane recovery (LPG product)'
    def _init(self, split=None):
        if split is None:
            split = {
                'Propane': 0.999,   # to LPG product (outs[0])
                'Hydrogen': 0.0,
                'CO2': 0.0,
                'CO': 0.0,
                'Water': 0.0,
                'Heptane': 0.0,
            }
        super()._init(split=split)
    def _setup(self):
        super()._setup()
        for s in self.outs:
            s.phase = 'g'  # or 'l' after HX302 if you merge chill into this block
    def _design(self):
        # Recovered propane in product stream (outs[0])
        self.design_results['Propane flow rate'] = self.outs[0].imass['Propane'] / 3600
           
            
            
@cost(
    'Total volume',ID='Protein extraction stage 1', units='m^3',cost=12080,S=1,n=0.525, CE=525.4,   BM=2.3,
   )
    
  #TODO update the cost and flow rate   
class ProteinExtractor1(bst.MixTank):
    """
    Stage-1 saline protein extraction mix tank (refer [27]).
    Meal + 0.1 M NaCl brine --> slurry to solids centrifuge (U203).
    Capital: BioSTEAM mix-tank correlation [Apostolakou 2009 / Seider].
    """
    line = 'Protein extraction stage 1'
    _default_tau = 2.0
    _default_kW_per_m3 = 1.0   #  TODO refer mixer_data agitator_kW_per_m3

    def _cost(self):
        self._decorated_cost()
        V = self.design_results.get('Total volume', 0)
        if V > 0:
            self.add_power_utility(self.kW_per_m3 * V)


@cost(
    'Total volume', ID='Protein extraction stage 2 wash', units='m^3', cost=12080, S=1,n=0.525, CE=525.4,BM=2.3,
)
  
  #TODO update the cost nand flow rate 
    
    
class ProteinExtractor2(bst.MixTank):
    """
    Stage-2 water wash mix tank (refer [27]).
    Pellet + wash water-->slurry to solids centrifuge (U205).
    """
    line = 'Protein extraction stage 2 wash'
    _default_tau = 1.0
    _default_kW_per_m3 = 1.0 #TODO refer mixer_data agitator_kW_per_m3

    def _cost(self):
        self._decorated_cost()
        V = self.design_results.get('Total volume', 0)
        if V > 0:
            self.add_power_utility(self.kW_per_m3 * V)            
            
            
# Total cost of wastewater treatment is combined into this placeholder
@cost(basis='Flow rate', ID='Wastewater system', units='kg/hr', 
      kW=7018.90125, S=393100, cost=50280080, CE=CEPCI[2010], n=0.6, BM=1)
class WastewaterSystemCost(Unit): pass       
        
        
            
            
            
            
            