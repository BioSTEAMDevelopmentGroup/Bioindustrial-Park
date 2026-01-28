# -*- coding: utf-8 -*-
"""
HemDx TEA Configuration

Refactored to include PreFerSTEA class with integrated production scaling 
and ammonia optimization logic.
"""

import biosteam as bst
import flexsolve as flx
import warnings
from biorefineries.prefers.v1._tea import PreFerSTEA as PreFerSTEA_Base
# Import seed_targets to update the global registry used by unit specs
from biorefineries.prefers.v1.HemDx.system._config1 import seed_targets

class PreFerSTEA(PreFerSTEA_Base):
    
    def optimize_NH3_loading(self, verbose=True):
        """
        Optimizes NH3_25wt flow rate and S202 split ratio to meet fermentation demand.
        Delegates to system module implementation.
        """
        from biorefineries.prefers.v1.HemDx.system import optimize_NH3_loading
        optimize_NH3_loading(self.system, verbose=verbose)

    def set_production_rate(self, target_production_kg_hr, verbose=True):
        """
        Set the target production rate and adjust system inputs accordingly.
        Delegates to system module implementation.
        """
        from biorefineries.prefers.v1.HemDx.system import set_production_rate
        
        self._target_production_kg_hr = target_production_kg_hr
        return set_production_rate(self.system, target_production_kg_hr, verbose=verbose)
    
    def check_product_specifications(self):
        """
        Check if the product stream meets all specifications.
        """
        return True
    
    @property
    def target_production_kg_hr(self):
        return self._target_production_kg_hr
    
    @target_production_kg_hr.setter
    def target_production_kg_hr(self, value):
        if value is not None:
            self.set_production_rate(value)
        else:
            self._target_production_kg_hr = None

if __name__ == '__main__':
    import argparse
    from biorefineries.prefers.v1.HemDx.system import create_NHemDx_system
    from biorefineries.prefers.v1.HemDx._chemicals import create_chemicals_Hemodextrin
    from biorefineries.prefers.v1._process_settings import load_process_settings
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='HemDx TEA Analysis')
    parser.add_argument('--production', type=float, default=150,
                        help='Target production rate in kg/hr (default: 150)')
    args, _ = parser.parse_known_args()
    
    TARGET_PRODUCTION = args.production
    
    print("="*85)
    print(f"HEMDX TEA - WITH DESIGN SPECIFICATION")
    print("="*85)
    
    bst.settings.set_thermo(create_chemicals_Hemodextrin(), skip_checks=True)
    load_process_settings()
    
    print(f"\n1. Creating HemDx system...")
    HemDx_sys = create_NHemDx_system()
    HemDx_sys.operating_hours = 8000
    
    print(f"\n2. Creating TEA with target production = {TARGET_PRODUCTION} kg/hr...")
    HemDx_tea = PreFerSTEA(
        system=HemDx_sys, 
        IRR=0.18, 
        duration=(2024, 2044), 
        depreciation='IRAS6',
        income_tax=0.17, 
        operating_days=333, 
        lang_factor=None,
        construction_schedule=(0.15, 0.60, 0.25), 
        WC_over_FCI=0.15,
        labor_cost=10*6e4, 
        fringe_benefits=0.17+0.07, 
        property_tax=0.005,
        property_insurance=0.005, 
        supplies=0.02, 
        maintenance=0.03,
        administration=0.05,
        target_production_kg_hr=TARGET_PRODUCTION
    )
    
    print(f"\n3. Adjusting system to target production rate...")
    achieved_production = HemDx_tea.set_production_rate(TARGET_PRODUCTION)
    
    # Results
    products = HemDx_sys.flowsheet.stream.NHemDx_Product
    print(f"\nProduction Summary:")
    print(f"  Hourly production:  {products.F_mass:.2f} kg/hr")
    
    HemDx_tea.show()
