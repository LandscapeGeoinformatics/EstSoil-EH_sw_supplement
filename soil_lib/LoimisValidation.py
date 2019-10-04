#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
/***************************************************************************

 Contains lookup functions for soil texture analysis
                             -------------------
        copyright            : (C) 2018-2019 by Alexander Kmoch
        email                : alexander.kmoch at ut.ee
 ***************************************************************************/
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the MIT License                                 *
 *                                                                         *
 ***************************************************************************/
"""

import os
import re
import sys
import traceback

import datetime
import logging

import math

from typing import (Any, Callable, Dict, Generic, Iterable, List, Mapping,
                    NewType, Sequence, Tuple, TypeVar, Union)

from operator import itemgetter
from subprocess import PIPE, Popen, call

import csv
import json

import numpy as np  # type: ignore
import pandas as pd  # type: ignore

import fiona  # type: ignore
import geopandas as gpd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

from .LoimisLookups import siffer_rules_lookup, error_lookup, texture_rules, kores_amp_strength, swat_ext_defaults_lookup
from .LoimisLookups import default_SOL_BD, default_SOL_CBN, default_SOL_AWC, default_SOL_K, default_SOL_ALB, default_USLE_K, default_ANION_EXCL, default_SOL_EC, default_SOL_CRK

from shapely.geometry import Point
from fiona.crs import from_epsg

from shapely import speedups

#################
# create logger
#################
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

logfile = logging.FileHandler(
    'temp_out/soil_convert.log', "w", encoding="utf-8")
console = logging.StreamHandler(sys.stdout)

# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logfile.setFormatter(formatter)
console.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(console)
logger.addHandler(logfile)

# #############################
# global vars and path names
##############################

speedups.enable()

soil_legend_file = "../soil/eesti_soil/mullalegend_4.csv"
mulla_df = pd.read_csv(
    soil_legend_file, encoding='latin1', sep=';', quotechar='"')
mulla_lookup = mulla_df['Tähistus kaardil'].tolist()
legend = set(mulla_lookup)


# Mahuprotsente ei kasutata korese tüüp 1, 9, 10, 11, ja 12 puhul
coarse = {'kr': {1: 'kruus'},
          'r':  {2: 'rähk'},
          'v':  {3: 'paeveeris'},
          'v0': {4: 'raudkiviveeris'},
          'kb': {5: 'klibu'},
          'ck': {6: 'kiltkivirähk'},
          'k':  {7: 'paekivid'},
          'k0': {8: 'raudkivid'},
          'pk': {9: 'paeplaadid'},
          'p':  {10: 'paas'},
          'd':  {11: 'liivakivi'},
          'lu': {12: 'lubisetted'}
          }


def joined_to_file():

    # read_excel Alar / encoding='latin1'
    emu_alar_df = pd.read_excel(
        r'C:\dev\05_geodata\soil\eesti_soil\alar_tartu_county_soil_decoded2016.xlsx')

    fields_to_check = ['MAPID', 'ID1', 'BONITEET', 'SIFFER', 'tüüpA', 'LOIMIS1',
                       'LA1', 'LA1_alg_cm', 'LA1_lopp_cm', 'korLA1', 'k_asteLA1a',
                       'LA2', 'LA2_alg_cm', 'LA2_lõpp_cm', 'korLA2', 'k_asteLA2a',
                       'LA3', 'LA3_alg_cm', 'LA3_lõpp_cm', 'korLA3', 'k_asteLA3a', 'X', 'Y']

    emu_alar_df = emu_alar_df[fields_to_check]
    logger.info(emu_alar_df.head(20))

    # tüüpA is our simplified_siffer (+/- in the end of tüüpA might be removed because of filtering in excel only)
    # M1 M2 M3 = M' M'' or M''' supposedly, but also AM777 been seen
    # LA1 LA2 is main texture loimis, but with __ for subscripts

    # load validation shape
    textures_df = gpd.read_file(
        'temp_out/eesti_soil_red1_texture_overview.shp')
    # textures_df = gpd.read_file(r'C:\dev\05_geodata\soil\eesti_soil\Mullakaart.shp')
    # textures_df = gpd.read_file(r'C:\dev\05_geodata\soil\eesti_soil\mapinfo_version\Mullakaart.TAB')

    fields = ['orig_fid', 'simplified', 'Boniteet', 'Loimis1', 'parse_info',
              'NLAYERS', 'SOL_ZMX',
              'EST_TXT1', 'SOL_Z1',
              'EST_TXT2', 'SOL_Z2',
              'EST_TXT3', 'SOL_Z3', 'geometry']

    textures_df = textures_df[fields]
    logger.info(textures_df.head(20))

    emu_alar_df['geometry'] = emu_alar_df.apply(
        lambda x: Point(x['X'], x['Y']), axis=1)
    geo_alar = gpd.GeoDataFrame(
        emu_alar_df, crs=from_epsg(3301), geometry='geometry')
    textures_df_aligned = textures_df.to_crs(geo_alar.crs)
    del(textures_df)

    # joined_df.to_file(driver='ESRI Shapefile', filename='temp_out/joined_df.shp', encoding="utf-8")
    # geo_alar.to_file(driver='ESRI Shapefile', filename='temp_out/alar_centroids_joined.shp', encoding="utf-8")

    joined_df = gpd.sjoin(textures_df_aligned, geo_alar,
                          how="inner", op="contains")
    logger.info(joined_df.index.size)
    logger.info(joined_df.dtypes)
    logger.info(joined_df.head())

    joined_df.to_file(driver='ESRI Shapefile',
                      filename='temp_out/joined_texture_validation_df.shp', encoding="utf-8")


# stufff
def replace_last(source_string, replace_what, replace_with):
    head, _sep, tail = source_string.rpartition(replace_what)
    return head + replace_with + tail


def adjust_topmarks(value):
    try:
        str_value = str(value)

        if str_value.endswith("'''"):
            return replace_last(str_value, "'''", '3')
        elif str(str_value).endswith("''"):
            return replace_last(str_value, "''", '2')
        elif str(str_value).endswith("'"):
            return replace_last(str_value, "'", '1')
        else:
            return str_value
    except TypeError as te:
        logger.error(value)
        logger.error(te)
        return value


def adjust_plusminus(row):
    if str(row['tüüpA']).endswith("-"):
        return replace_last(row['tüüpA'], '-', '')
    elif str(row['tüüpA']).endswith("+"):
        return replace_last(row['tüüpA'], '+', '')
    else:
        return row['tüüpA']


def comp_sifs(x):
    a = str(x['simplified_compare'])
    b = str(x['tyypA_compare'])
    if a == b:
        return 1
    else:
        return 0


def startswithshort(row):
    rests = ['TxR', 'TzM', 'TxM', 'M', 'AM', "S", 'LkIg']
    retval = 0
    for t in rests:
        if str(row['tüüpA']).startswith(t) or str(row['tüüpA']) == "" or str(row['tüüpA']) == "nan":
            retval = 1
    return retval


def adjust_texture_subs(x):
    return x.replace("₁", "_1").replace("₂", "_2").replace("₃", "_3").replace("⁰", "_")


def comp_txts(x):
    a = str(x['parse_trans'])
    b = str(x['liiv_loimis'])
    if (pd.isna(a) or a == "no_peenes" or a is None) and (pd.isna(b) or b == "nan" or b is None):
        return 1
    if a == b:
        return 1
    else:
        return 0


def comp_txt1(x):
    a = x['EST_TXT1']
    b = x['LA1']
    if (pd.isna(a) or a == "no_peenes" or a is None) and (pd.isna(b) or b == "nan" or b is None):
        return 1
    if str(a) == str(b):
        return 1
    else:
        return 0


def comp_txt2(x):
    a = x['EST_TXT2']
    b = x['LA2']
    if (pd.isna(a) or a == "no_peenes" or a is None or a == "None") and (pd.isna(b) or b == "nan" or b is None):
        return 1
    if a == b:
        return 1
    else:
        return 0


def comp_txt3(x):
    a = x['EST_TXT3']
    b = x['LA3']
    if (pd.isna(a) or a == "no_peenes" or a is None) and (pd.isna(b) or b == "nan" or b is None):
        return 1
    if a == b:
        return 1
    else:
        return 0


def compare_soil_types():
    joined_df = gpd.read_file('temp_out/joined_texture_validation_df.shp')

    logger.info(joined_df.dtypes)
    logger.info('## sol types siffer')

    joined_df['simplified_compare'] = joined_df['simplified'].apply(
        lambda x: adjust_topmarks(x))
    joined_df['tyypA_compare'] = joined_df.apply(
        lambda x: adjust_plusminus(x), axis=1)
    joined_df['comp_sifs'] = joined_df.apply(lambda x: comp_sifs(x), axis=1)

    logger.info("## joined_df['comp_sifs'].sum()")
    logger.info(joined_df['comp_sifs'].sum())

    neg_df = joined_df.loc[joined_df['comp_sifs'] == 0].copy()
    logger.info('## neg_df.index.size')
    logger.info(neg_df.index.size)
    logger.info('## neg_df.sample(2)')
    logger.info(neg_df.sample(2))
    neg_df['rests'] = neg_df.apply(lambda x: startswithshort(x), axis=1)
    neg_df = neg_df[neg_df['rests'] != 1]

    logger.info('## neg_df.index.size')
    logger.info(neg_df.index.size)
    logger.info("## neg_df['comp_sifs'].sum()")
    logger.info(neg_df['comp_sifs'].sum())

    summary_df1 = neg_df.groupby(['simplified_compare', 'tyypA_compare']).size(
    ).reset_index().rename(columns={0: 'count'})
    logger.info('## ummary_df1.sample(2)')
    logger.info(summary_df1.sample(2))
    logger.info('## summary_df1.index.size')
    logger.info(summary_df1.index.size)
    summary_df1.to_csv(
        'temp_out/soil_compare_siffer_summary_df.csv', encoding='utf-8')

    logger.info('## update_main_siffer, find_main_siffer')

    from SWAT_paragen_dask import update_main_siffer, find_main_siffer
    demo_short_sif = neg_df[['orig_fid',
                             'simplified_compare', 'tyypA_compare']].copy()
    demo_short_sif.rename(columns={'tyypA_compare': 'Siffer'}, inplace=True)
    demo_short_sif['Sif1'] = demo_short_sif['Siffer']
    demo_short_sif[['tmp_fix_siffer', 'siffer_match_info']] = demo_short_sif.apply(
        lambda x: update_main_siffer(x), axis=1)

    logger.info('## demo_short_sif.sample(2) with tmp_fix_siffer')
    logger.info(demo_short_sif.sample(2))

    demo_short_sif = demo_short_sif[demo_short_sif['siffer_match_info']
                                    == "not_matched"]
    logger.info("## demo_short_sif.index.size not_matched")
    logger.info(demo_short_sif.index.size)
    logger.info('## demo_short_sif. with tmp_fix_siffer and not_matched')
    logger.info(demo_short_sif)

    demo_short_sif = demo_short_sif[demo_short_sif['siffer_match_info']
                                    == "none_error"]
    logger.info("## demo_short_sif.index.size none_error")
    logger.info(demo_short_sif.index.size)

    demo_short_sif.rename(
        columns={'Siffer': 'simplified_compare', 'Sif1': 'tyypA_compare'}, inplace=True)
    summary_df1_s = demo_short_sif.groupby(['tmp_fix_siffer', 'tyypA_compare']).size(
    ).reset_index().rename(columns={0: 'count'})
    logger.info("## summary_df1_s after siffer match")
    logger.info(summary_df1_s)
    logger.info("## summary_df1_s.index.size")
    logger.info(summary_df1_s.index.size)
    summary_df1_s.to_csv(
        'temp_out/soil_compare_siffer_double_simplified_summary_df.csv', encoding='utf-8')
    joined_df.to_file(driver='ESRI Shapefile',
                      filename='temp_out/joined_texture_validation_siffers_df.shp', encoding="utf-8")


def compare_textures():

    joined_df = gpd.read_file(
        'temp_out/joined_texture_validation_siffers_df.shp')

    logger.info(joined_df.dtypes)
    logger.info('## texture')
    # texture

    joined_df['parse_trans'] = joined_df['parse_info'].apply(
        lambda x: adjust_texture_subs(str(x)))
    joined_df['liiv_loimis'] = joined_df['LOIMIS1_1'].apply(
        lambda x: str(x).replace("liiv", "l").replace("savi", "s"))

    joined_df['comp_txts'] = joined_df.apply(lambda x: comp_txts(x), axis=1)

    logger.info("## joined_df.sample(2)")
    joined_df[['Loimis1', 'parse_info', 'parse_trans',
               'liiv_loimis', 'ID1']].sample(2)
    logger.info("## joined_df['comp_txts'].sum()")
    logger.info(joined_df['comp_txts'].sum())

    logger.info("## comp_txt1 comp_txt2 comp_txt3")
    joined_df['comp_txt1'] = joined_df.apply(lambda x: comp_txt1(x), axis=1)
    joined_df['comp_txt2'] = joined_df.apply(lambda x: comp_txt2(x), axis=1)
    joined_df['comp_txt3'] = joined_df.apply(lambda x: comp_txt3(x), axis=1)

    logger.info("## joined_df['comp_txt1'].sum()")
    logger.info(joined_df['comp_txt1'].sum())

    logger.info(joined_df.loc[joined_df['comp_txt1'] == 0].fillna('no_value').groupby(
        ['EST_TXT1',  'LA1']).size().reset_index().rename(columns={0: 'count'}))
    logger.info(joined_df.loc[joined_df['comp_txt1'] == 0].fillna('no_value').groupby(
        ['EST_TXT1',  'LA1']).size().reset_index().rename(columns={0: 'count'})['count'].sum())
    joined_df.loc[joined_df['comp_txt1'] == 0].fillna('no_value').groupby(['EST_TXT1',  'LA1']).size().reset_index(
    ).rename(columns={0: 'count'}).to_csv('temp_out/soil_compare_texture_summary_layer_1.csv', encoding='utf-8')

    logger.info("## joined_df['comp_txt2'].sum()")
    logger.info(joined_df['comp_txt2'].sum())

    logger.info(joined_df.loc[joined_df['comp_txt2'] == 0].fillna('no_value').groupby(
        ['EST_TXT2',  'LA2']).size().reset_index().rename(columns={0: 'count'}))
    logger.info(joined_df.loc[joined_df['comp_txt2'] == 0].fillna('no_value').groupby(
        ['EST_TXT2',  'LA2']).size().reset_index().rename(columns={0: 'count'})['count'].sum())
    joined_df.loc[joined_df['comp_txt2'] == 0].fillna('no_value').groupby(['EST_TXT2',  'LA2']).size().reset_index(
    ).rename(columns={0: 'count'}).to_csv('temp_out/soil_compare_texture_summary_layer_2.csv', encoding='utf-8')

    logger.info("## joined_df['comp_txt3'].sum()")
    logger.info(joined_df['comp_txt3'].sum())
    logger.info(joined_df.loc[joined_df['comp_txt3'] == 0].fillna('no_value').groupby(
        ['EST_TXT3',  'LA3']).size().reset_index().rename(columns={0: 'count'}))
    logger.info(joined_df.loc[joined_df['comp_txt3'] == 0].fillna('no_value').groupby(
        ['EST_TXT3',  'LA3']).size().reset_index().rename(columns={0: 'count'})['count'].sum())
    joined_df.loc[joined_df['comp_txt3'] == 0].fillna('no_value').groupby(['EST_TXT3',  'LA3']).size().reset_index(
    ).rename(columns={0: 'count'}).to_csv('temp_out/soil_compare_texture_summary_layer_3.csv', encoding='utf-8')

    logger.info("## joined_df['comp_txt ... '].sum() without None filled")
    logger.info(joined_df.loc[joined_df['comp_txt1'] == 0].groupby(
        ['EST_TXT1',  'LA1']).size().sum())
    logger.info(joined_df.loc[joined_df['comp_txt2'] == 0].groupby(
        ['EST_TXT2',  'LA2']).size().sum())
    logger.info(joined_df.loc[joined_df['comp_txt3'] == 0].groupby(
        ['EST_TXT3',  'LA3']).size().sum())

    joined_df.to_file(driver='ESRI Shapefile',
                      filename='temp_out/joined_texture_validation_textures_df.shp', encoding="utf-8")


def write_out_textures_for_ksat_rosetta():
    
    ksat_rosetta_src = gpd.read_file(
        'temp_out/eesti_soil_red1_swat_ext_values.shp')

    logger.info(ksat_rosetta_src.dtypes)
    logger.info('## Ksat Rosetta loader')

    # 2 NEW SSC (sand, silt, clay)
    fields_min = ['SOL_SAND1', 'SOL_SILT1', 'SOL_CLAY1',
                  'SOL_SAND2', 'SOL_SILT2', 'SOL_CLAY2',
                  'SOL_SAND3', 'SOL_SILT3', 'SOL_CLAY3']

    # 3 NEW SSC BD (sand, silt, clay, bulk density)
    fields_bd = ['SOL_SAND1', 'SOL_SILT1', 'SOL_CLAY1', 'SOL_BD1',
                 'SOL_SAND2', 'SOL_SILT2', 'SOL_CLAY2', 'SOL_BD2',
                 'SOL_SAND3', 'SOL_SILT3', 'SOL_CLAY3', 'SOL_BD3']

    overall_rows = ksat_rosetta_src.index.size
    logger.info(overall_rows)

    # UNITS!
    # SSC in weight %
    # Rosetta ANN needs BD in g/cm3 -> we provide between 1.1 and 1.9 Mg/m³ (seems to fit M mega)
    ## we don't have TH33 and T1500 as cm3/cm3

    for step in [100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000]:

        logger.info("{} {}".format(fields_min[0:3], step))
        x1 = ksat_rosetta_src.as_matrix(columns=fields_min[0:3])
        np.savetxt("temp_out/Rosetta-3.0beta/output/mulla_ksat_lay1_{}_mod2_input_np.txt".format(step), x1[step-100000:step], delimiter="  ")

        logger.info("{} {}".format(fields_min[3:6], step))
        x1 = ksat_rosetta_src.as_matrix(columns=fields_min[3:6])
        np.savetxt("temp_out/Rosetta-3.0beta/output/mulla_ksat_lay2_{}_mod2_input_np.txt".format(step), x1[step-100000:step], delimiter="  ")

        logger.info("{} {}".format(fields_min[6:9], step))
        x1 = ksat_rosetta_src.as_matrix(columns=fields_min[6:9])
        np.savetxt("temp_out/Rosetta-3.0beta/output/mulla_ksat_lay3_{}_mod2_input_np.txt".format(step), x1[step-100000:step], delimiter="  ")

        logger.info("{} {}".format(fields_min[0:4], step))
        x2 = ksat_rosetta_src.as_matrix(columns=fields_bd[0:4])
        np.savetxt( "temp_out/Rosetta-3.0beta/output/mulla_ksat_lay1_{}_mod3_input_np.txt".format(step), x2[step-100000:step], delimiter="  ")

        logger.info("{} {}".format(fields_min[4:8], step))
        x2 = ksat_rosetta_src.as_matrix(columns=fields_bd[4:8])
        np.savetxt( "temp_out/Rosetta-3.0beta/output/mulla_ksat_lay2_{}_mod3_input_np.txt".format(step), x2[step-100000:step], delimiter="  ")

        logger.info("{} {}".format(fields_min[8:12], step))
        x2 = ksat_rosetta_src.as_matrix(columns=fields_bd[8:12])
        np.savetxt( "temp_out/Rosetta-3.0beta/output/mulla_ksat_lay3_{}_mod3_input_np.txt".format(step), x2[step-100000:step], delimiter="  ")


def yield_series_for_layers():

    # we'll get 6 series, K1, K2, K3 (Ksat per layer), as model 2 (SSC) and model 3 (SSC BD)

    # OUTPUT
    # theta_r [cm3/cm3]
    # theta_s [cm3/cm3]
    # alpha  [1/cm]
    # n
    # Ks in [cm/day] -> SWAT needs mm/hr
    # standard deviations apply to the log10 forms for alpha, n and KS
    # NOT their their antilog forms 

    for layer in range(1,4):

        x_mod2_all = np.zeros(0,)
        x_mod3_all = np.zeros(0,)
        for step in [100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000]:
        
            logger.info("mod2 layer {} step {}".format(layer, step))
            x2 = np.loadtxt("temp_out/Rosetta-3.0beta/output/output_mulla_ksat_lay{}_{}_mod2_input_np.txt".format(layer,step), delimiter=",")
            col_ksat = x2[:,4]
            x_mod2_all = np.append(x_mod2_all, col_ksat)

            logger.info("mod3 layer {} step {}".format(layer, step))
            x3 = np.loadtxt( "temp_out/Rosetta-3.0beta/output/output_mulla_ksat_lay{}_{}_mod3_input_np.txt".format(layer,step), delimiter=",")
            col_ksat = x3[:,4]
            x_mod3_all = np.append(x_mod3_all, col_ksat)
        
        # apply div 2.4 for cm/day -> mm/hr
        yield (pd.Series(x_mod2_all / 2.4), pd.Series(x_mod3_all / 2.4))


def read_in_ksat_predicted():
    # we'll get 6 series, K1, K2, K3 (Ksat per layer), as model 2 (SSC) and model 3 (SSC BD)
    ksat_rosetta_repl = gpd.read_file(
        'temp_out/eesti_soil_red1_swat_ext_values.shp')
    
    layer = 1
    for series_tup in yield_series_for_layers():
        ksat_ssc, ksat_ssc_bd = series_tup

        ksat_field = 'SOL_K{}_new'.format(layer)
        ksat_field_bd = 'SOL_K{}_new_BD'.format(layer)
        ksat_rosetta_repl[ksat_field] = ksat_ssc
        ksat_rosetta_repl[ksat_field_bd] = ksat_ssc_bd
        layer = layer + 1
        if layer >= 4:
            break
    
    ksat_rosetta_repl.rename(columns={'SOL_K1': 'SOL_K1_old','SOL_K2': 'SOL_K2_old','SOL_K3': 'SOL_K3_old'}, inplace=True)

    csv_drop1 = ksat_rosetta_repl[['orig_fid',
        'SOL_K1_old','SOL_K1_new','SOL_K1_new_BD',
        'SOL_K2_old','SOL_K2_new','SOL_K2_new_BD',
        'SOL_K3_old','SOL_K3_new','SOL_K3_new_BD' ]]
    csv_drop1.to_csv('temp_out/eesti_soil_red1_swat_ext_values_rosetta_ksat.csv', encoding="utf-8")
    del(csv_drop1)

    ksat_rosetta_repl.rename(columns={'SOL_K1_new': 'SOL_K1','SOL_K2_new': 'SOL_K2','SOL_K3_new': 'SOL_K3'}, inplace=True)
    ksat_rosetta_repl.drop(columns=['SOL_K1_old','SOL_K1_new_BD',
        'SOL_K2_old','SOL_K2_new_BD',
        'SOL_K3_old','SOL_K3_new_BD' ], inplace=True)
    ksat_rosetta_repl.to_file(driver='ESRI Shapefile', filename='temp_out/eesti_soil_red1_swat_ext_values_new_ksat.shp', encoding="utf-8")


if __name__ == '__main__':

    plt.style.use('ggplot')
    plt.rcParams['figure.figsize'] = (25, 10)
    pd.set_option('display.max_rows', 100)
    pd.set_option('display.max_columns', 100)

    mp.freeze_support()

    num_cores: int = int(mp.cpu_count()-2)
    start = datetime.datetime.now()  # strftime("%Y-%m-%d %H:%M:%S")
    global_time_counter = start
    logger.info(
        "initialising at ... {} with num cores: {}".format(start.strftime("%Y-%m-%d %H:%M:%S"), num_cores))

    # joined_to_file()
    # compare_soil_types()
    # compare_textures()
    # write_out_textures_for_ksat_rosetta()
    ## Rpredict.py
    ## model 2 SSC predict_for_mod2.bat
    ## model 3 SSC BD predict_for_mod3.bat
    read_in_ksat_predicted()

    # Soil grids
    #
    # load and put in shape summary
    #
    # https://gis.stackexchange.com/questions/260304/extract-raster-values-within-shapefile-with-pygeoprocessing-or-gdal/260380
    #
    # C:\dev\05_geodata\soil\soilgrids_download
