# EstSoil-EH_sw_supplement

[![DOI](https://zenodo.org/badge/212613251.svg)](https://zenodo.org/badge/latestdoi/212613251)

For more information on the development of this dataset look for "EstSoil-EH: a high-resolution eco-hydrological modelling parameters dataset for Estonia", Alexander Kmoch, Arno Kanal†, Alar Astover, Ain Kull, Holger Virro, Aveliina Helm, Meelis Pärtel, Ivika Ostonen and Evelyn Uuemaa, 2021, Earth Syst. Sci. Data, 13, 83–97, [https://doi.org/10.5194/essd-13-83-2021](https://doi.org/10.5194/essd-13-83-2021)

Supplementary Jupyter notebooks, external scripts, and GIS application usages and codes


- steps to mostly seamlessly derive most values via pandas/geopandas dataframe functions ('soil_lib' python module, '01_soilmap_soiltypes_textures_layers.ipynb')
  - the standardising of soil types plus WRB soil types
  - the grammar, parsing and standadising loimis/texture codes
  - deriving layers, depths, sand, silt, clay and rock content
- external step to predict K_sat with the Rosetta 3 program ('Rosetta-3.0') 
- additional step to load the predicted K_sat values from the Rosetta 3 program (05_hydrogrids_extents_and_awc_extract.ipynb)
- external steps to calculated slope, LS, TRI and TWI based on Estonian 5m DEM and summarized (mean) into soil units with QGIS 3.4
- external step to create SOC soil samples training dataset ('03_SOC_RF_preps.ipynb')
- additional steps to derive soil prarameters
  - loading the SOC sample data, hyperparameter tuning and applying RF model to predict SOC for all soil units ('04_soilmap_SOC-RF_BD.ipynb')
  - calculating BD based on SOC ('04_soilmap_SOC-RF_BD.ipynb')
  - additional step to calculate AWC raster from FC and WP rasters from EU-HydroSoilGrids v1.0 (https://zenodo.org/record/3446747)
    and summarise into soil units for correct depths ('05_hydrogrids_extents_and_awc_extract.ipynb')
  - additional step to calculate USLE_K ('04_soilmap_SOC-RF_BD.ipynb')

- pip/conda requirements/environment yml file for estsoil grammar and lookups library and archived jupyter notebooks for reference: requirements.txt
- Rosetta 3 program: only worked on Python 2.7, conda environment requirements.txt included in Rosetta-3.0 folder/zip

