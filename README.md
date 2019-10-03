# EstSoil-EH_sw_supplement

Supplementary Jupyter notebooks, external scripts, and GIS application usages and codes


- steps to mostly seamlessly derive most values via pandas/geopandas dataframe functions (soil_lib python module, soilmap_swat_gen_fasttrack.ipynb)
  - the standardising of soil types plus WRB soil types
  - the grammar, parsing and standadising loimis/texture codes
  - deriving layers, depths, sand, silt, clay and rock content
- external step to predict K_sat with the Rosetta 3 program
- additional step to load the predicted K_sat values from the Rosetta 3 program ("soil update Ksat and AWC work.ipynb")
- external steps to calculated slope, LS, TRI and TWI based on Estonian 5m DEM and summarized (mean) into soil units with QGIS 3.4
- external step to create SOC soil samples training dataset (soilmap_soc_rf_preps.ipynb, SOC-RF.ipynb)
- additional steps to derive soil prarameters
  - loading the SOC sample data, hyperparameter tuning and applying RF model to predict SOC for all soil units (SOC-RF.ipynb)
  - calculating BD based on SOC (SOC-RF.ipynb)
  - additional step to calculate AWC raster from FC and WP rasters from EU-HydroSoilGrids v1.0 (https://zenodo.org/record/3446747)
    and summarise into soil units for correct depths (rasterstats_hydrogrids_redo_awc.ipynb)
  - additional step to calculate USLE_K ("soil update Ksat and AWC work.ipynb")

- pip/conda requirements/environment yml file for estsoil grammar and lookups library and archived jupyter notebooks for reference: requirements.txt
- Rosetta 3 program: only worked on Python 2.7, conda environment requirements.txt included in Rosetta-3.0 folder/zip

