# ResinColumnAdv Specification

**Advanced Resin Column Unit** with mechanistic modeling and bidirectional parameter modes.

Supported by module: `biorefineries.prefers.v1._units_adv`

## Overview
`ResinColumnAdv` is a standalone unit operation built on `bst.Unit` that unifies Ion Exchange and Adsorption modeling. It extends the original `ResinColumn` functionality with:
1.  **Mechanistic Performance Models**: Yield and impurity removal are functions of column volumes (CVs) via exponential decay curves.
2.  **Efficiency Factor**: A global scalar (`efficiency`, default 0.85) to represent real-world deviations from ideal binding/elution.
3.  **Bidirectional Parameter Modes**: 
    - `'ratio'` mode: Users set performance targets (e.g., Yield=0.99) directly.
    - `'CV'` mode: Users set operating parameters (CVs), and performance is calculated dynamically.

## Parameters

### Modes
- `preset`: 'IonExchange' or 'Adsorption'
- `parameter_mode`: 'ratio' (default) or 'CV'

### Physical Parameters
| Parameter                     | Default | Description                                                                                   |
| :---------------------------- | :------ | :-------------------------------------------------------------------------------------------- |
| `operating_pressure_drop_bar` | 2.0     | Pressure drop across the column (used for pump sizing).                                       |
| `efficiency`                  | 0.85    | Efficiency scalar (0.0-1.0). In 'ratio' mode, usually set to 1.0 if ratios are final targets. |

### Performance Parameters (Ratio Mode Targets)
| Parameter                  | Description                                                             |
| :------------------------- | :---------------------------------------------------------------------- |
| `TargetProduct_Yield`      | Fraction of target recovered in Eluate/Product stream.                  |
| `BoundImpurity_Removal`    | (IEX Only) Fraction of bound impurities removed to Regeneration stream. |
| `NonBinding_Carryover`     | (IEX Only) Fraction of non-binding species entrained in Product stream. |
| `Wash_Impurity_Carryover`  | (Adsorption Only) Fraction of impurities removed to Wash stream.        |
| `Regen_Impurity_Carryover` | (Adsorption Only) Fraction of impurities removed to Regen stream.       |

### Operating Parameters (CV Mode Inputs)
| Parameter         | Description                                                          |
| :---------------- | :------------------------------------------------------------------- |
| `elution_CV`      | Column volumes of elution buffer. Increases Yield.                   |
| `wash_CV`         | Column volumes of wash buffer. Increases Impurity Removal (to Wash). |
| `regeneration_CV` | Column volumes of regen buffer. Increases Recovery/removal to Regen. |

## Mechanistic Models (CV Mode)

### Adsorption Preset
- **Yield**: $Yield = Yield_{base} \times \eta \times (1 - e^{-k_{elute} \cdot CV_{elute}})$
- **Impurity Removal (Wash)**: $R_{wash} = R_{base,wash} \times \eta \times (1 - e^{-k_{wash} \cdot CV_{wash}})$
- **Impurity Removal (Regen)**: $R_{regen} = R_{base,regen} \times \eta \times (1 - e^{-k_{regen} \cdot CV_{regen}})$

Default calibration assumes industrial standard performance at typical CVs (e.g., 99% recovery at ~3 CV).

## Streams
1.  **ins[0]**: Feed
2.  **ins[1]**: Wash Buffer
3.  **ins[2]**: Elution Buffer
4.  **ins[3]**: Regeneration Buffer

1.  **outs[0]**: Flowthrough / Treated Waste
2.  **outs[1]**: Eluate / Product
3.  **outs[2]**: Regeneration Waste
4.  **outs[3]**: Spent Wash Buffer
