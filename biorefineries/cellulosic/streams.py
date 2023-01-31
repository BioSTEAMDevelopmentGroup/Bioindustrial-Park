# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""
from biosteam import stream_kwargs

cornstover = stream_kwargs(
    'cornstover',
    Glucan=0.28,
    Xylan=0.1562,
    Galactan=0.001144,
    Arabinan=0.01904,
    Mannan=0.0048,
    Lignin=0.12608,
    Acetate=0.01448,
    Protein=0.0248,
    Extract=0.1172,
    Ash=0.03944,
    Sucrose=0.00616,
    Water=0.2,
    total_flow=104229.16,
    units='kg/hr',
    price=0.05158816935126135,
)

switchgrass = stream_kwargs(
    'switchgrass',
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
    units='kg/hr',
    price=0.08, # Price of switchgrass, table 4, Madhu Khanna et al. (Costs of producing miscanthus and switchgrass for bioenergy in Illinois);  https://www.sciencedirect.com/science/article/pii/S096195340700205X?casa_token=KfYfzJtDwv0AAAAA:OqeJmpofk1kIgFk2DcUvXNG35qYwlWvPKZ7ENI3R6RUKeoahiTDpOhhd_mpLtRthTGuXJKDzMOc
)

NaOH = stream_kwargs(
    'NaOH', 
    NaOH=1, 
    price=0.14952852932689323
)

ammonia = stream_kwargs(
    'ammonia', 
    phase='l',
    NH3=1, 
    P=12 * 101325, 
    price=0.4485966110937889,
)

sulfuric_acid = stream_kwargs(
    'sulfuric_acid',
    P=5.4*101325,
    T=294.15,
    Water=130,
    H2SO4=1800,
    units='kg/hr',
    price=0.0897171175961359
)
denaturant = stream_kwargs(
    'denaturant',
    Octane=1,
    price=0.756,
)
ethanol = stream_kwargs(
    'ethanol',
    price=0.7198608114634679,
)
lignin = stream_kwargs('lignin')
slurry = stream_kwargs('slurry')
DAP = stream_kwargs('DAP', price=0.9869213628968231)
CSL = stream_kwargs('CSL', price=0.0568241480781522)
vent = stream_kwargs('vent')
beer = stream_kwargs('beer')
pretreated_biomass = stream_kwargs('pretreated_biomass')
cellulase = stream_kwargs('cellulase', price=0.212)
saccharification_water = stream_kwargs('saccharification_water')
slurry = stream_kwargs('slurry')
nanofilter_retentate = stream_kwargs('nanofilter_retentate')
pretreatment_wastewater = stream_kwargs('pretreatment_wastewater')
juice = stream_kwargs('juice')
hydrolysate = stream_kwargs('hydrolysate')