## Reference data

This directory contains the reference data as share in Ranc et al. 2022. The code base will be tested against these values, i.e. the code should approximately provide the same results given the same drivers, parameters and reference values.

### Drivers

Driver data consists of:

- slope (EU DEM + quadratic term; slope.asc, slope_sq.asc)
- tree cover (at a 325m resolution multi-grain SSA - corine land cover + quadratic term, tcd_325grain.asc, tcd_325grain_sq.asc)
- land cover (reforested; agriculture and pastures maps at 0.05 ha, landcover_5322.asc and landcover_agri.asc)

### Tracks

Roe deer tracks for the full set of observations (`Aspromonte_roedeer_traj.txt`) and a subset for a single individual 1196 (`Aspromonte_roedeer_traj_1196.txt`).

### Config

The configuration files of different runs of the original publication, only the `config_best_Mmem_fitting.txt` should be used as a reference set of parameters.

### Validation

The output of model runs (`objective_function.*`), and a composite `global_resources.asc` which is the composite of the driver data and determines the order of the layers in the data cube (generated from the driver  maps).

