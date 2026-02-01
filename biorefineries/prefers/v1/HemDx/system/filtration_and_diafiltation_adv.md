# Mechanistic Models for Filtration & Diafiltration (TEA/LCA Focus)

## Overview
This document outlines "Backward-Calculation" mechanistic models for filtration unit operations. These models are designed for Techno-Economic Analysis (TEA) and Life Cycle Assessment (LCA) to bridge the gap between process targets (Yield, Purity) and engineering requirements (Membrane Area, Buffer Usage, Power).

---

## 1. Tangential Flow Filtration (TFF) / Diafiltration
**Purpose:** Washing out impurities (buffer exchange) or concentrating product.
**Operation Mode:** Continuous loop (feed -> membrane -> retentate recycle).

### A. Diafiltration (Buffer Exchange) Model
This model calculates the required buffer volume to achieve a specific impurity removal target. It assumes a CSTR (Continuous Stirred-Tank Reactor) washout mechanism.

**Backward Formula (Target $\rightarrow$ Design):**
To calculate the number of Diavolumes ($N_{DV}$) required:

$$N_{DV} = \frac{\ln \left( \frac{C_{imp, initial}}{C_{imp, final}} \right)}{1 - \sigma_{imp}}$$

* **$C_{imp, final} / C_{imp, initial}$**: The target fraction of impurity remaining (e.g., 0.01 for 99% removal).
* **$\sigma_{imp}$ (Rejection Coefficient)**: 
    * $\sigma = 0$: Freely passing (Salts, water).
    * $\sigma = 1$: Fully retained (Product).
    * *Typical for salts:* 0.0 – 0.05.

**Buffer Volume Calculation:**
$$V_{buffer} = N_{DV} \times V_{current\_retentate}$$

### B. Concentration (Ultrafiltration) Model
This model calculates the volume reduction required to hit a concentration target.

**Backward Formula:**
$$VCF = \left( \frac{C_{target, final}}{C_{target, initial}} \right)^{1/\sigma_{target}}$$

* **$VCF$ (Volumetric Concentration Factor)**: $V_{feed} / V_{retentate}$.
* **$\sigma_{target}$**: Usually assumed 1.0 (100% retention) for the product.

**Permeate Volume to Remove:**
$$V_{permeate} = V_{feed} \times \left( 1 - \frac{1}{VCF} \right)$$

---

## 2. The Physics of Flux (Sizing the Hardware)
To calculate **Membrane Area**, you need the average Flux ($J$). In TEA, simple constant flux is often inaccurate for high-concentration processes.

### The Gel Polarization Model (Limiting Flux)
As concentration increases, flux decreases due to the formation of a gel layer on the membrane wall.

$$J(C) = k \cdot \ln\left( \frac{C_{gel}}{C_{bulk}} \right)$$

* **$C_{gel}$**: The concentration of the product at the membrane wall (saturation limit).
* **$k$**: Mass transfer coefficient (dependent on pump crossflow velocity).

### Average Flux Approximation
For the backward calculation, integrate the flux over the concentration step:

$$J_{avg} \approx J_{initial} - k \left( 1 - \frac{1}{VCF} \right)$$

### Hardware Sizing Formula
$$Area_{m2} = \frac{V_{total\_permeate} (L)}{J_{avg} (LMH) \times t_{process} (h)}$$

---

## 3. Dead-End Filtration (Clarification)
**Purpose:** Removing solids/cells.
**Operation Mode:** Single pass, limit is capacity (clogging).

**Backward Formula:**
$$Area_{needed} = \frac{V_{batch}}{V_{max} \times \eta_{safety}}$$

* **$V_{max}$ (L/m²)**: Maximum capacity before filter plugs (experimental parameter).
* **$\eta_{safety}$**: Safety margin (typically 0.7 - 0.8).

---

## 4. Python Implementation (Ready for BioSTEAM)

This Python class performs the backward calculation: input your targets, output your hardware and utility requirements.

```python
import numpy as np

class TFF_Unit_Sizer:
    def __init__(self, design_flux_LMH, rejection_imp=0.1, rejection_target=1.0):
        """
        design_flux_LMH: Expected average flux (L/m2/h)
        rejection_imp: 0.0 = free pass, 1.0 = full retention
        """
        self.J = design_flux_LMH
        self.sigma_imp = rejection_imp
        self.sigma_target = rejection_target

    def backward_design(self, vol_feed_L, target_purity_removal, target_concentration_factor):
        """
        Reverse calculations for TFF sizing.
        
        Args:
            vol_feed_L: Initial batch volume
            target_purity_removal: Fraction of impurity to remove (e.g. 0.99 for 99%)
            target_concentration_factor: VCF (e.g. 10x concentration)
            
        Returns:
            Dictionary containing hardware size (Area) and Utility (Buffer)
        """
        
        # --- Step 1: Concentration (Volume Reduction) ---
        # We concentrate first to save buffer in the next step
        vol_retentate_L = vol_feed_L / target_concentration_factor
        vol_permeate_conc_L = vol_feed_L - vol_retentate_L
        
        # --- Step 2: Diafiltration (Washing) ---
        # Calculate Diavolumes (N) needed to hit purity target
        # Formula: N = ln(C_in / C_out) / (1 - sigma)
        
        fraction_remaining = 1.0 - target_purity_removal
        
        if fraction_remaining <= 0:
             raise ValueError("Purity removal must be < 1.0")
             
        # Guard against division by zero if sigma is 1
        denom = 1.0 - self.sigma_imp
        if denom < 1e-6:
             N_dv = 0 # Cannot wash out if fully retained
        else:
             N_dv = np.log(1.0 / fraction_remaining) / denom
        
        vol_buffer_L = N_dv * vol_retentate_L
        
        # --- Step 3: Sizing (Membrane Area) ---
        # Total liquid volume pushed through the membrane
        total_vol_permeate = vol_permeate_conc_L + vol_buffer_L
        
        # Constraint: Process must finish in X hours (e.g., 4h shift)
        process_time_h = 4.0 
        area_m2 = total_vol_permeate / (self.J * process_time_h)
        
        return {
            "Membrane_Area_m2": area_m2,
            "Buffer_Usage_L": vol_buffer_L,
            "Total_Permeate_Waste_L": total_vol_permeate,
            "Diavolumes": N_dv,
            "Final_Retentate_Vol_L": vol_retentate_L
        }

# --- Example Usage ---
if __name__ == "__main__":
    # Setup: 50 LMH flux, Salts pass freely (0.05 rejection)
    sizer = TFF_Unit_Sizer(design_flux_LMH=50, rejection_imp=0.05)
    
    # Scenario: 1000L Feed, Want 99% Salt Removal, 5x Concentration
    results = sizer.backward_design(
        vol_feed_L=1000, 
        target_purity_removal=0.99, 
        target_concentration_factor=5.0
    )

    print("--- Backward Design Results ---")
    print(f"Required Membrane Area: {results['Membrane_Area_m2']:.2f} m²")
    print(f"Buffer Required:        {results['Buffer_Usage_L']:.2f} L")
    print(f"Diavolumes Needed:      {results['Diavolumes']:.2f} DV")