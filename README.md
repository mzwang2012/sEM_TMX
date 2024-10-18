# sEM_TMX
Stochastic emulator of daily maximum temperature

Reference:
Stochastic Emulators of Spatially Resolved Extreme Temperatures of Earth System Models
by Mengze Wang, Andre Souza, Raffaele Ferrari, and Themistoklis Sapsis
https://doi.org/10.22541/essoar.172858084.46299070/v1


This emulator is constructed using CMIP6 data. Two CMIP6 models have been attemped: (1) CNRM-CM6-1-HR model; (2) MPI-ESM1-2-LR model
These data can be retrieved through the Earth System Grid Federation interface https://aims2.llnl.gov/search/cmip6/.

1. File structure
   The folder structure and files under CNRM/ and MPI/ are almost the same. The difference in the code occurs because (1) The grid resolution of CNRM and MPI models are different; (2) MPI has a large ensemble of datasets while CNRM only has one realization; (3) We investigated the influence of ensemble size for MPI model but not for CNRM.

2. Under either CNRM/ or MPI/, first run the codes under EOF/ to compute the Empirical Orthogonal Functions.

    2.1 Run nc_global_tas_mean.m to compute the global mean temperature

    2.2 Run pca_global_historical.m to get the EOFs, using historical data.

    2.3 Run pca_global_coeff_historical.m to project the TMX data onto different EOFs to get the time series of EOF coefficients.

    2.4 Run pca_global_post.m to plot the leading EOFs and the total variance represented by first N modes.

    2.5 Run pca_global_coeff_post.m to plot the time series of EOF coefficients

3. Under Emulator/

   3.1 First run stats_global_coarse.m to generate the true statistics, including the 10-year summer mean, standard deviation, and 97.5% quantile. These stats are stored to disk.

   3.2 Run stats_global_corase_rom_01.m to generate emulated statistics. When this code is run, it calls rom_glbT_lin_mean_var_pred.m to generate emulated time series of EOF coefficients. These coefficients are then combined with EOFs to get instantaneous TMX and then compute the emulated statistics.
