# -*- coding: utf-8 -*-
"""
/***************************************************************************

 Summarizes raster statistics in parallel from soilgrids to estonian soil polygons
                             -------------------
        copyright            : (C) 2018-2020 by Alexander Kmoch
        email                : alexander.kmoch at ut.ee
 ***************************************************************************/
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the MIT License                                 *
 *                                                                         *
 ***************************************************************************/
"""
from dask.distributed import Client
import geopandas as gpd

import pandas as pd

import fiona

from rasterstats import zonal_stats

import logging
import datetime

log_level = logging.INFO
# create logger
logger = logging.getLogger(__name__)
logger.setLevel(log_level)

fh = logging.FileHandler('script_output.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)

start = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

if __name__ == "__main__":
    # client = Client()
    client = Client(processes=True,
                    n_workers=2,
                    threads_per_worker=1,
                    memory_limit='12GB')

    print(client.scheduler_info()['services'])
    logger.info("client ready at ... {} ... at {}".format(client.scheduler_info()['services'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    soilgrids = {
        'sand': 'soil250_grid_sand_sd',
        'silt': 'soil250_grid_silt_sd',
        'clay': 'soil250_grid_clay_sd',
        'rock': 'soil250_grid_coarsefrag_sd',
        'bd': 'soil250_grid_bulkdens_sd',
        'soc': 'soil250_grid_soc_sd',
        'awc': 'soil250_grid_awc_sd',
        'k_sat': 'soil250_grid_k_sat_sd'
    }

    raster_file_collection = []

    template_raster_conf = {
        'variable_name': 'sand',
        'layer_num': 1,
        'actual_file_ref': '/run/user/1817077476/alex_tmp_geo/soilgrids_download/soil250_grid_sand_sd6_3301.tif'
        # 'actual_file_ref': '../soilgrids_download/soil250_grid_sand_sd6_3301.tif'
    }

    base_path = '/home/DOMENIS.UT.EE/kmoch/soil_paper_materials/WORK_TMP'
    for layer_num in range(1,8):
    
        for layer_type in soilgrids.keys():
            file_name = f"{base_path}/soilgrids_download/{soilgrids[layer_type]}{layer_num}_3301.tif"
            # file_name = f"../soilgrids_download/{soilgrids[layer_type]}{layer_num}_3301.tif"
            try:
                with open(file_name, 'r') as fh:
                    # Load configuration file values
                    print(file_name)
                    template_raster_conf = {
                        'variable_name': layer_type,
                        'layer_num': layer_num,
                        'actual_file_ref': file_name
                    }
                    raster_file_collection.append(template_raster_conf)
            except FileNotFoundError:
                # Keep preset values
                print("NOT " + f"{soilgrids[layer_type]}{layer_num}")
                logger.warn("NOT " + f"{soilgrids[layer_type]}{layer_num}")

    # '../data_deposit/EstSoil-EH_id_tmp.shp'
    geom_template = f'{base_path}/EstSoil-EH_id_tmp.shp'
    
    print("reading in template ESTSOIL geoms")
    logger.info("reading in template ESTSOIL geoms ... {} ... at {}".format(geom_template, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    
    geom_template_df = gpd.read_file(geom_template, encoding='utf-8')

    refdf_scattered = client.scatter(geom_template_df, broadcast=True)

    # for raster_conf_dict in raster_file_collection:
        
    def inner_raster_summary(raster_conf_dict, scattered_df):
    
        # raster_summary(raster_conf_dict)
        # print(f"Starting with {raster_conf_dict['actual_file_ref']}")
        print("Starting with {} at {}".format(raster_conf_dict['actual_file_ref'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        logger.info("Starting with {} at {}".format(raster_conf_dict['actual_file_ref'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        
        variable_name = raster_conf_dict['variable_name']
        layer_num = raster_conf_dict['layer_num']
        actual_file_ref = raster_conf_dict['actual_file_ref']
            
        tif_src = actual_file_ref
            
        # shp_tmp_out = f"../data_deposit/EstSoil-EH_{variable_name}{layer_num}_zonal_stats.shp"
        csv_tmp_out = f"../data_deposit/EstSoil-EH_{variable_name}{layer_num}_zonal_stats.csv"
        parquet_tmp_out = f"../data_deposit/EstSoil-EH_{variable_name}{layer_num}_zonal_stats.parquet.gzip"
        
        with fiona.open(geom_template) as vector_src:
        
            # src_crs = vector_src.crs
            # display(src_crs)
            # src_schema = vector_src.schema
            # display(src_schema)
            # src_schema['properties']["mean"] = "float"
            # src_schema['properties']["std"] = "float"
            
            print("zonal stats")
            outputs = zonal_stats(vector_src,
                    tif_src,
                    stats="mean std",
                    # geojson_out=True,
                    all_touched=True)
            
            # with fiona.open(shp_tmp_out, "w", driver="ESRI Shapefile", schema=src_schema, crs=src_crs) as collection:
            #     collection.writerecords(outputs)
            #     print(len(collection))
            #     collection.flush()
            
            print("output df preps")
            geo_stats_df = pd.DataFrame(outputs)
            print(geo_stats_df.columns)
            
            var1_name = f"TEST_{variable_name.upper()}{layer_num}"
            var2_name = f"TEST_{variable_name.upper()}{layer_num}_STD"
            
            geo_stats_df.rename(columns={'mean': var1_name,'std': var2_name}, inplace=True)
        
            print(geo_stats_df.columns)

            geom_df = pd.concat([scattered_df, geo_stats_df], axis=1)
            geom_df.drop(columns=['geometry'], inplace=True)
            geom_df = geom_df[['orig_fid', var1_name, var2_name]]
            geom_df.set_index('orig_fid', inplace=True)
        
            print("writing {} at {}".format(parquet_tmp_out, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            logger.info("writing {} at {}".format(parquet_tmp_out, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            
            geom_df.to_csv(csv_tmp_out, encoding='utf-8')
            geom_df.to_parquet(parquet_tmp_out, compression='gzip', engine='pyarrow')
            
        del(geom_df)
        del(outputs)
        del(geo_stats_df)
            
        return f"{csv_tmp_out} written"

    futures = [client.submit(inner_raster_summary, raster_conf_obj, refdf_scattered) for raster_conf_obj in raster_file_collection]
    # futures = client.map(inner_raster_summary, raster_file_collection)(refdf_scattered)
    results = client.gather(futures)
    for i in results:
        print(i)

    # for raster_conf_dict in raster_file_collection:



