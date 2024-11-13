# User Instructions
1. For installation and spike-sorting, please read the manual.
2. For basic analysis, combining "Presentation" log files and matlab-based unit files, use codes in the folder analysis_postphy_core.
3. The first code to run is "pool_data_ephys", which combines all selected unit files into a single Matlab dictionary. Currently the code is set to pool data into recorded regions.
4. The next code to run is "gather_raster_ephys", which will convert spike-times into rasters and rates.
5. The code "plot_raster" is optional, but will allow plotting of individual units as a spike-raster and lick-raster plot.   
