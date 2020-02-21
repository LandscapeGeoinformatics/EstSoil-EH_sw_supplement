#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
/***************************************************************************

 parallel loimis test, fix, and parse
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


from typing import (Any, Callable, Dict, Generic, Iterable, List, Mapping,
                    NewType, Sequence, Tuple, TypeVar, Union)

import os
import sys
import traceback
import logging
import datetime

import dask.dataframe as dd
from dask.distributed import Client

import pandas as pd

from arpeggio import ParserPython, PTNodeVisitor, visit_parse_tree, NoMatch
from arpeggio import RegExMatch, Optional, ZeroOrMore, OneOrMore, EOF, UnorderedGroup, And, Not, Combine


from soil_lib.LoimisLookups import siffer_rules_lookup, updated_texture_error_lookup, fillers_by_numbers,texture_rules

from soil_lib.LoimisGrammarV2 import (update_main_siffer_lt, split_and_cut_dask_sharp, new_grammar,
                                      test_brackets_dask, parse_test_dask_multiple,
                                      parse_reconstituate_dask, load_default_texture_defensively_dask )
                                      
from soil_lib.LoimisVisitor import (LpVisitor, can_parse_multiple_get_parser, loimis_grammar_product_dask_multiple,
                                    test_layer_depths_dask, set_texture_values_dask)



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

logger.info("starting now {}".format(start))
print("starting now {}".format(start))

parquet_loimis_unique = "loimis_unique.parquet.gzip"

if __name__ == "__main__":
    # client = Client("127.0.0.1:63288")
    # client = Client(processes=True,
    #                n_workers=6,
    #                 threads_per_worker=1,
    #                 memory_limit='1GB')

    # print(client.scheduler_info()['services'])
    # logger.info("client ready at ... {} ... at {}".format(
    #         client.scheduler_info()['services'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    
    # ddf = dd.read_sql_table("landboard_soilmap_orig_2017_fixedgeom",
    #                     PG_DB_URL,
    #                     "orig_fid", columns=['orig_fid','Loimis1','Varv'],
    #                     schema="geoworkspace" )
    ddf = dd.read_parquet(parquet_loimis_unique, engine='pyarrow')

    ddf = ddf.repartition(npartitions=30)
    
    result = ddf.map_partitions(
    lambda pdf: pdf['Loimis1'].apply(
        lambda x: split_and_cut_dask_sharp(x, updated_texture_error_lookup)),
        meta=('split_layered', 'str')
    )
    
    ddf = ddf.assign(split_layered=result)
    
    
    result = ddf.map_partitions(
    lambda pdf: pdf['split_layered'].apply(
        lambda x: test_brackets_dask(x)),
        meta=('num_brackets_fixed', 'str' )
    )
    
    ddf = ddf.assign(num_brackets_fixed=result)
    
    
    def per_partition_lambda_parse_test(p_df):
        parsers_d = new_grammar()
        return p_df['num_brackets_fixed'].apply(lambda x: parse_test_dask_multiple(x, parsers_d))
    
    
    result = ddf.map_partitions(per_partition_lambda_parse_test, meta=pd.DataFrame({0: ['str'], 1: [ 1 ] }))
    
    result = result.rename(columns={0: 'test_parse', 1: 'test_parse_errors'})
    
    ddf = ddf.assign(test_parse=result['test_parse'], test_parse_errors=result['test_parse_errors'])
    
    
    def per_partition_lambda_parse_reconstituate(p_df):
        return p_df['test_parse'].apply(lambda x: parse_reconstituate_dask(x))
    
    
    result = ddf.map_partitions(per_partition_lambda_parse_reconstituate, meta=pd.DataFrame({0: ['str'], 1: [ 1 ] }))
    
    result = result.rename(columns={0: 'loimis_reconst', 1: 'has_no_info'})
    
    ddf = ddf.assign(loimis_reconst=result['loimis_reconst'], has_no_info=result['has_no_info'])
    
    ddf = ddf.persist()
    
    
    def per_partition_lambda_loimis_grammar(p_df):
        parsers_d = new_grammar()
        return p_df['loimis_reconst'].apply(lambda x:  loimis_grammar_product_dask_multiple(x, parsers_d))


    result = ddf.map_partitions(per_partition_lambda_loimis_grammar, meta=pd.DataFrame({0: [{'a':'b'}], 1: ['str'] }))
    
    result = result.rename(columns={0: 'loimis_grammar', 1: 'parse_info'})
    
    ddf = ddf.assign(loimis_grammar=result['loimis_grammar'], parse_info=result['parse_info'])


    def per_partition_lambda_layer_depths(p_df):
        return p_df['loimis_grammar'].apply(lambda x: test_layer_depths_dask(x))
    
    
    data_struct_meta = { 0: [0],  1: [0],  2: [0],  3: [0],  4: [0],  5: [0] }
    
    result = ddf.map_partitions(per_partition_lambda_layer_depths, meta=pd.DataFrame(data_struct_meta))
    
    cols = { 0: 'nlayers',  1: 'SOL_ZMX',  2: 'SOL_Z1',  3: 'SOL_Z2',  4: 'SOL_Z3',  5: 'SOL_Z4' }
    
    result = result.rename(columns=cols)
    
    ddf = ddf.assign(nlayers=result['nlayers'], SOL_ZMX=result['SOL_ZMX'],
                    SOL_Z1=result['SOL_Z1'], SOL_Z2=result['SOL_Z2'],
                    SOL_Z3=result['SOL_Z3'], SOL_Z4=result['SOL_Z4'])


    def per_partition_lambda_texture_values(p_df):
        return p_df['loimis_grammar'].apply(lambda x: set_texture_values_dask(x))
    
    
    data_struct_meta = { 0: ['str'],  1: ['str'],  2: ['str'],  3: [0],  4: [0],  5: [0],  6: [0],
                         7: ['str'],  8: ['str'],  9: ['str'], 10: [0], 11: [0], 12: [0], 13: [0],
                        14: ['str'], 15: ['str'], 16: ['str'], 17: [0], 18: [0], 19: [0], 20: [0],
                        21: ['str'], 22: ['str'], 23: ['str'], 24: [0], 25: [0], 26: [0], 27: [0] , 28: [{'a': 'str'}]}
    
    result = ddf.map_partitions(per_partition_lambda_texture_values,
                                meta=pd.DataFrame(data_struct_meta))
    
    cols = { 0: 'EST_TXT1',  1: 'EST_CRS1',  2: 'LXTYPE1',  3: 'SOL_CLAY1',  4: 'SOL_SILT1',  5: 'SOL_SAND1',  6: 'SOL_ROCK1',
             7: 'EST_TXT2',  8: 'EST_CRS2',  9: 'LXTYPE2', 10: 'SOL_CLAY2', 11: 'SOL_SILT2', 12: 'SOL_SAND2', 13: 'SOL_ROCK2',
            14: 'EST_TXT3', 15: 'EST_CRS3', 16: 'LXTYPE3', 17: 'SOL_CLAY3', 18: 'SOL_SILT3', 19: 'SOL_SAND3', 20: 'SOL_ROCK3',
            21: 'EST_TXT4', 22: 'EST_CRS4', 23: 'LXTYPE4', 24: 'SOL_CLAY4', 25: 'SOL_SILT4', 26: 'SOL_SAND4', 27: 'SOL_ROCK4', 28: 'loimis_search' }
    
    result = result.rename(columns=cols)
    
    ddf = ddf.assign( EST_TXT1=result['EST_TXT1'], EST_CRS1=result['EST_CRS1'], LXTYPE1=result['LXTYPE1'], SOL_CLAY1=result['SOL_CLAY1'], SOL_SILT1=result['SOL_SILT1'], SOL_SAND1=result['SOL_SAND1'], SOL_ROCK1=result['SOL_ROCK1'],
                    EST_TXT2=result['EST_TXT2'], EST_CRS2=result['EST_CRS2'], LXTYPE2=result['LXTYPE2'], SOL_CLAY2=result['SOL_CLAY2'], SOL_SILT2=result['SOL_SILT2'], SOL_SAND2=result['SOL_SAND2'], SOL_ROCK2=result['SOL_ROCK2'],
                    EST_TXT3=result['EST_TXT3'], EST_CRS3=result['EST_CRS3'], LXTYPE3=result['LXTYPE3'], SOL_CLAY3=result['SOL_CLAY3'], SOL_SILT3=result['SOL_SILT3'], SOL_SAND3=result['SOL_SAND3'], SOL_ROCK3=result['SOL_ROCK3'],
                    EST_TXT4=result['EST_TXT4'], EST_CRS4=result['EST_CRS4'], LXTYPE4=result['LXTYPE4'], SOL_CLAY4=result['SOL_CLAY4'], SOL_SILT4=result['SOL_SILT4'], SOL_SAND4=result['SOL_SAND4'], SOL_ROCK4=result['SOL_ROCK4'], loimis_search=result['loimis_search'] )
    

    
    