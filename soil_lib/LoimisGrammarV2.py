# -*- coding: utf-8 -*-
"""
/***************************************************************************

 Contains text analytical functions for soil texture analysis
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

from typing import (Any, Callable, Dict, Generic, Iterable, List, Mapping,
                    NewType, Sequence, Tuple, TypeVar, Union)

import pandas as pd
import re

import logging

from arpeggio import (ParserPython, EOF, And, Combine, NoMatch, Optional,
                      RegExMatch, ZeroOrMore, Not, OneOrMore)

from .LoimisLookups import siffer_rules_lookup, updated_texture_error_lookup, fillers_by_numbers

# create logger
logger = logging.getLogger(__name__)


##############################
# Arpeggio-based Grammar for the Texture value encoding (loimis)
#################################

def new_grammar():
    
    # "+" r, v, kb,k, pk, p, lu - karbonaatne(kihisev materjal), but seems also to be used in more contexts
    def kPlus(): return '+'
    
    # can only be one digit (should be subscript, otherwise error), diff to depth numbers
    # def amp1(): return Combine(['₁', (Not(depth_range), '1')])
    
    # def amp2(): return Combine(['₂', (Not(depth_range), '2')])
    
    # def amp3(): return Combine(['₃', (Not(depth_range), '3')])
    
    def amp1(): return '₁'
    
    def amp2(): return '₂'
    
    def amp3(): return '₃'
    
    def amp4(): return '₄'
    
    def amp5(): return '₅'
    
    def amplifiers(): return [amp1, amp2, amp3, amp4, amp5]
    
    
    def depth_number(): return RegExMatch(r'\d+')  # at least two digits for depth numbers (well might as well be one)
    
    def depth_range(): return Optional(kPlus), depth_number, ZeroOrMore('-', depth_number)
    
    
    # peenes
    def l(): return 'l'  # liiv, en: sand
    
    def pl(): return 'pl'  # pl - peenliiv, en: fine sand (täiendina peenliivakas)
    
    def plsl(): return 'plsl'  # plsl - peenliivakas saviliiv, en: fine clayey sand
    
    def sl(): return 'sl'  # sl - saviliiv, en: clayey sand
    
    def tsl(): return 'tsl'  # tsl - tolmjas saviliiv, en: dusty clayey sand
    
    def dk(): return 'dk'  # dk - liivakivirähk, en: sandstone grit
    
    def ls(): return 'ls'
    
    def tls(): return 'tls'  # tls - tolmjas liivsavi, en: dusty sandy clay
    
    def s(): return 's'  # s - savi, en: clay
    
    def peenes_list(): return [ dk, sl, ls, s, l ]
    
    def peenes_list_t(): return [ tsl, tls ]
    
    def peenes_list_p(): return [ plsl, pl ]
    
    def peenes(): return Optional(kPlus), peenes_list, Optional(amplifiers), Optional(depth_range)
    
    def peenes_t(): return Optional(kPlus), peenes_list_t, Optional(amplifiers), Optional(depth_range)
    
    def peenes_p(): return Optional(kPlus), peenes_list_p, Optional(amplifiers), Optional(depth_range)
    
    def th(): return 'th'  # th 15 või th 15-20 toorhuumusliku horisondi tüsedus, en: raw humus thickness
    
    # def t(): return Combine( [ (Not(tls), 't'), (Not(tsl), 't')] )
    def t(): return 't'
    
    def th_turfs(): return Optional(kPlus), th, Optional(amplifiers), Optional(depth_range)
    
    def t_turfs(): return Optional(kPlus), t, Optional(amplifiers), Optional(depth_range)
    
    # kores
    # kr, p, d, lu, pk - astmeid ei kasutata
    
    def kr(): return 'kr'  # kr - kruus, en: gravel
    
    def p(): return 'p' # p - paas, en: limestone rock (p_l_ does not exist)
    
    def d(): return 'd'  # d - liivakivi, en: sandstone rock?
    
    def lu(): return 'lu'  # lu - lubisetted, calcic deposits? kores peenes turfs alternateComma amplifiers amp1 depth_range
    
    def pk(): return 'pk'  # pk - paeplaadid, en: limestone rock plates?
    
    def ck(): return 'ck'  # ck - kiltkivirähk, en: debris, rubble of slate/schist
    
    def no_info(): return 'no_info'
    
    def skeleton_no_amp(): return [  pk, kr, d, lu, ck, no_info ]
    
    # r, v, kb, k (1-5, or none)
    def r_0(): return ['r⁰', 'r°']  # r - rähk, en: grit, rubble
    
    def r(): return 'r'  # r - rähk, en: grit, rubble r⁰
    
    def r_norm(): return [ r_0, r ]
    
    def v_0(): return ['v⁰', 'v°']  # v° - raudkiviveeris, en: granite rock
    
    def v(): return 'v'  # v - paeveeris, en: limestone rock
    
    def v_norm(): return [ v_0, v]
    
    def kb(): return 'kb'  # kb - klibu, en: debris, rubble
    
    def kb_0(): return ['kb⁰', 'kb°']
    
    def kb_norm(): return [ kb_0, kb ]
    
    def k_0(): return ['k⁰', 'k°']  # k⁰ - raudkivid, en: granite rock
    
    def k(): return 'k'  # k - paekivid, en: limestone rock
    
    def k_norm(): return [ k_0, k ]
    
    def skeleton_for_amp(): return [  r_norm, v_norm, kb_norm, k_norm, no_info ]
    
    
    def amp_kores(): return Optional(kPlus), skeleton_for_amp, Optional(amplifiers), Optional(depth_range)
    
    def no_amp_kores(): return Optional(kPlus), skeleton_no_amp, Optional(depth_range)
    
    def p_kores(): return Optional(kPlus), p, Optional(depth_range)
    
    def double_kores_consitutent_1(): return ZeroOrMore(no_amp_kores), Optional(th_turfs), [peenes_p, peenes_t]
    def double_kores_consitutent_2(): return ZeroOrMore(no_amp_kores), Optional(t_turfs), Optional(peenes)
    def double_kores_consitutent_3(): return ZeroOrMore(no_amp_kores), Optional(t_turfs), peenes_p
    def double_kores_consitutent_4(): return ZeroOrMore(no_amp_kores), Optional(th_turfs), Optional(peenes)
    
    def double_kores_consitutent_1a(): return ZeroOrMore(amp_kores), Optional(th_turfs), [peenes_p, peenes_t]
    def double_kores_consitutent_2a(): return ZeroOrMore(amp_kores), Optional(t_turfs), Optional(peenes)
    def double_kores_consitutent_3a(): return ZeroOrMore(amp_kores), Optional(t_turfs), peenes_p
    def double_kores_consitutent_4a(): return ZeroOrMore(amp_kores), Optional(th_turfs), Optional(peenes)
    
    def double_kores_consitutent_5(): return ZeroOrMore(p_kores), Optional(th_turfs), [peenes, peenes_t]
    def double_kores_consitutent_6(): return ZeroOrMore(p_kores), Optional(t_turfs), Optional(peenes)
    
    
    def loimisp_short_dk_1(): return double_kores_consitutent_1, EOF
    def loimisp_short_dk_2(): return double_kores_consitutent_2, EOF
    def loimisp_short_dk_3(): return double_kores_consitutent_3, EOF
    def loimisp_short_dk_4(): return double_kores_consitutent_4, EOF
    def loimisp_short_dk_1a(): return double_kores_consitutent_1a, EOF
    def loimisp_short_dk_2a(): return double_kores_consitutent_2a, EOF
    def loimisp_short_dk_3a(): return double_kores_consitutent_3a, EOF
    def loimisp_short_dk_4a(): return double_kores_consitutent_4a, EOF
    def loimisp_short_dk_5(): return double_kores_consitutent_5, EOF
    def loimisp_short_dk_6(): return double_kores_consitutent_6, EOF
    
    sp1 = ParserPython(loimisp_short_dk_1,memoization=False)
    sp2 = ParserPython(loimisp_short_dk_2,memoization=False)
    sp3 = ParserPython(loimisp_short_dk_3,memoization=False)
    sp4 = ParserPython(loimisp_short_dk_4,memoization=False)
    sp1a = ParserPython(loimisp_short_dk_1a,memoization=False)
    sp2a = ParserPython(loimisp_short_dk_2a,memoization=False)
    sp3a = ParserPython(loimisp_short_dk_3a,memoization=False)
    sp4a = ParserPython(loimisp_short_dk_4a,memoization=False)
    sp5 = ParserPython(loimisp_short_dk_5,memoization=False)
    sp6 = ParserPython(loimisp_short_dk_6,memoization=False)
    
    parsers_d = {'sp1': sp1,
               'sp2':sp2,
               'sp3':sp3,
               'sp4':sp4,
               'sp1a':sp1a,
               'sp2a':sp2a,
               'sp3a':sp3a,
               'sp4a':sp4a,
               'sp5':sp5,
               'sp6':sp6 }
    
    return parsers_d

###################################
# commonly used regular expressions
###################################

bracket_matcher_inside = re.compile("\((?P<inside>.*)\)", re.IGNORECASE)

bracket_matcher = re.compile("\(.*\)", re.IGNORECASE)

bracket_matcher_num1 = re.compile("(?P<keep>\d+)\(\d+\)") # 90(40)

bracket_matcher_num2 = re.compile("\((?P<keep>\d+-\d+)\)") # (70-100)

bracket_matcher_num3 = re.compile("(?P<keep>\d+-\d+)\(\d+\)") # 30-50(100)

bracket_matcher_num4 = re.compile("(?P<keep>\d+-\d+)\(\d+-\d+\)") # 50-90(40-50)

bracket_matcher_num5 = re.compile("(?P<keep>\d+)\(\d+-\d+\)") # 90(40-50)

bracket_matcher_num6 = re.compile("\((?P<keep>\d+)\)") # (70)

bracket_matcher_num7 = re.compile("\(\d+\)(?P<keep>\d+-\d+)") # (40)50-90

open_bracket_matcher = re.compile("\(.*$", re.IGNORECASE)

subs_rep_matcher = re.compile("(?P<first>₁|₂|₃|₄|₅)(₁|₂|₃|₄|₅)*(\ |-|,|\+)(₁|₂|₃|₄|₅)")


##################################
# utility functions
#################################

def find_main_siffer(siffer_string: str, mulla_lookup: List) -> str:
    """finds the main soil type for a single soil type string

    Arguments:
        siffer_string {string} -- the soil type in Estonian with possible extended codes

    Returns:
        str -- the soil type in Estonian reduced to the allowed values from the legend
    """

    # always replace "’" with "'"
    siffer = siffer_string.replace("’", "'")
    foundFor = 0
    whereBefore = ''
    foundAny = False
    for sif_map in mulla_lookup:
        if siffer and isinstance(siffer, str) and siffer.startswith(sif_map) and len(sif_map) > foundFor:
            whereBefore = sif_map
            foundFor = len(sif_map)
            foundAny = True
    if foundAny == False:
        # logger.info('havent found any sif_map for %s' % (siffer))
        return 'not_matched'
    else:
        return whereBefore


def update_main_siffer_l(sif: str, mulla_lookup: List) -> str:

    # always replace "’" with "'"
    siffer = sif.replace("’", "'")
    result1 = find_main_siffer(siffer, mulla_lookup)
    if result1 == 'not_matched':
        try:
            siffer_lookup1 = siffer_rules_lookup[siffer]
            return pd.Series([siffer_lookup1, 'rules_lookup_siffer'])
        except KeyError as ex:
            return pd.Series([siffer, 'not_matched'])
    else:
        return pd.Series([result1, 'siffer'])


def update_main_siffer_lt(sif: str, mulla_lookup: List) -> pd.Series(Tuple):

    # always replace "’" with "'"
    siffer = sif.replace("’", "'")
    result1 = find_main_siffer(siffer, mulla_lookup)
    if result1 == 'not_matched':
        try:
            siffer_lookup1 = siffer_rules_lookup[siffer]
            return pd.Series((siffer_lookup1, 'rules_lookup_siffer'))
        except KeyError as ex:
            return pd.Series((siffer, 'not_matched'))
    else:
        return pd.Series((result1, 'siffer'))


def split_and_cut(l_str: str) -> List[str]:
    main_part_only: str = l_str.split(";")[0]
    try:
        main_part_only = updated_texture_error_lookup[main_part_only]
    except KeyError as ex:
        pass
    layered: List[str] = main_part_only.split("/")
    return layered


def split_and_cut_dask(l_str: str, updated_texture_error_lookup) -> Tuple[str]:
    main_part_only: str = l_str.split(";")[0]
    try:
        main_part_only = updated_texture_error_lookup[main_part_only]
    except KeyError as ex:
        pass
    layered: List[str] = main_part_only.split("/")
    return tuple(layered)


def split_and_cut_dask_sharp(l_str: str, updated_texture_error_lookup) -> str:
    main_part_only: str = l_str.split(";")[0]
    try:
        main_part_only = updated_texture_error_lookup[main_part_only]
    except KeyError as ex:
        pass
    layered: List[str] = main_part_only.split("/")
    return str('||'.join(layered))


def repl_bt_reg(rex: Any, txt: str) -> str:
    v = rex.search(txt)
    return txt.replace(v.group(), v.group('keep'))


def consolidate_num_bracket(e: str) -> str:
    # r1 90(40) & r3 30-50(100)
    # r2 (70-100) & r4 50-90(40-50)
    # r2 (70-100) & r5 90(40-50)
    # r3 50-90(40-50) & r5 90(40-50)
    
    # r3 before r5
    r3 = bracket_matcher_num3.search(e)
    if r3 is not None:
        return repl_bt_reg(bracket_matcher_num3, e)
    
    # r4 before r2
    r4 = bracket_matcher_num4.search(e)
    if r4 is not None:
        return repl_bt_reg(bracket_matcher_num4, e)
        
    # r5 before r2
    r5 = bracket_matcher_num5.search(e)
    if r5 is not None:
        return repl_bt_reg(bracket_matcher_num5, e)
    
    # r3 before r1
    r1 = bracket_matcher_num1.search(e)
    if r1 is not None:
        return repl_bt_reg(bracket_matcher_num1, e)
    
    r2 = bracket_matcher_num2.search(e)
    if r2 is not None:
        return repl_bt_reg(bracket_matcher_num2, e)
    
    r6 = bracket_matcher_num6.search(e)
    if r6 is not None:
        return repl_bt_reg(bracket_matcher_num6, e)
    
    return e


def subscripts_consolidate(e: str) -> str:
    txt = e
    if subs_rep_matcher.search(txt) is not None:
        v = subs_rep_matcher.findall(txt)
        # print(v)
        for g in v:
            rp = ''.join(g)
            txt = txt.replace(rp, g[0])

    return txt


def test_brackets(l_str: List[str]) -> List[str]:
    fixed_list = []
    for e in l_str:
        txt = subscripts_consolidate(e)
        result = bracket_matcher.search(txt)
        if result is not None:
            num_fix = consolidate_num_bracket(txt)
            # if num_fix == txt:
            #     global_brackets_list.append(txt)
            fixed_list.append(num_fix)
        else:
            fixed_list.append(txt)
    return fixed_list


def test_brackets_dask(l_str: str) -> str:
    fixed_list = []
    for e in l_str.split('||'):
        txt = subscripts_consolidate(e)
        result = bracket_matcher.search(txt)
        if result is not None:
            num_fix = consolidate_num_bracket(txt)
            # if num_fix == txt:
            #     global_brackets_list.append(txt)
            fixed_list.append(num_fix)
        else:
            fixed_list.append(txt)
    return "||".join(fixed_list)


def remove_brkt_inside(rex: Any, txt: str) -> str:
    v = rex.search(txt)
    return txt.replace(v.group(), v.group('inside'))


def remove_inside_group(rex: Any, txt: str) -> str:
    v = rex.search(txt)
    return txt.replace(v.group(), '')


def can_parse(e: str, parser: Any) -> bool:
    try:
        parse_tree = parser.parse(e)
        return True
    except NoMatch as nm:
        return False


def can_parse_multiple(e: str, parsers_d: Dict) -> bool:
    done = False
    for parser_name in parsers_d.keys():
        parser = parsers_d[parser_name]
        try:
            done = can_parse(e, parser)
        except Exception as ex:
            continue
        if done:
            break
    return done



def clean_dashes(elem: str):
    a1 = elem.replace('puudub', '').replace('<Null>', '').replace('tuhk', 'tls').replace('vesi', '')
    # .replace('Kog', '').
    a2 = a1.replace('\ufeff', '').replace('/0', '').replace('\x81', '').replace('la', '').replace('al', '').replace('üle', '+')
    
    a3 = a2 .replace('r-ls', 'rls').replace('ls/3','ls₃').replace('ls⁰', 'ls').replace('v⁰₁sl-v⁰₁ls','v⁰₁sl')
    a4 = a3.replace('sl-ls', 'sl').replace('tsl-tls', 'tsl').replace('v⁰₁l-v⁰₁sl', 'v⁰₁l').replace('l-sl', 'l').replace('sl-l','sl')
    
    # .replace('-ls', 'ls').replace('lls','ls').replace('lsl','ls₁').replace('lss','ls')
    
    a5 = a4.replace('₄₃', '₄')
    a6 = a5.replace('₁₁', '₁').replace('₁₂', '₁').replace('₂₁', '₂').replace('₁₃', '₁').replace('₃₁', '₃').replace('₃₅', '₃')
    a7 = a6.replace('₃₂', '₃').replace('₃₄', '₃').replace('₂₃', '₂').replace('₂₄', '₂').replace('₂₅', '₂').replace('₂₂', '₂')
    a8 = a7.replace('kr₁', 'kr').replace('kr₂', 'kr').replace('kr₃', 'kr').replace('kr₄', 'kr').replace('kr₅', 'kr') # .replace('_', '')
    
    
    a9 = a8.replace('kr₁', 'kr').replace('kr₂', 'kr').replace('kr₃', 'kr').replace('kr₄', 'kr').replace('kr₅', 'kr')
    
    a10 = a9.replace('ck₁', 'ck').replace('ck₂', 'kr').replace('ck₃', 'kr').replace('ck₄', 'kr').replace('ck₅', 'kr')
    a11 = a10.replace('e', '').replace('n', '').replace('i', '').replace('m', '')
    
    a12 = a11.replace('~', '').replace('%', '').replace('>', '+').replace('*', '').replace('++', '+')
    
    a13 = a12.replace('Ko', 'k⁰').replace('ko', 'k⁰').replace('o', '⁰').replace('K', '⁰').replace('⁰⁰','⁰').replace('⁰_f⁰', '')
    a14 = a13.replace('turvas', 't').replace('lubi', 'lu').replace('vee', '').replace('/mergel', '').replace('prügi', '').replace('killustik', 'ck')
    
    return a14


def strip_more(x: str) -> str:
    no_comma = x.split(",")
    bm = open_bracket_matcher.search(no_comma[0])
    nxt = ''
    if bm is not None:
        nxt = remove_inside_group(open_bracket_matcher, no_comma[0])
    else:
        nxt = no_comma[0]
    xdash = clean_dashes(nxt)
    return xdash


def consolidate_loimis(t_e: str, t_parser: Any, tk_parser: Any) -> (str, int):
    e = t_e
    try:
        e = updated_texture_error_lookup[e]
    except KeyError as ex:
        pass
    v = bracket_matcher_inside.search(e)
    # test if has brackets
    if v is not None:
        t1 = v.group('inside')
        # test if inside can be parsed at all
        if can_parse(t1, t_parser):
            # valid inside, test with inclusion in main string (just removing the brackets)
            v_1 = remove_brkt_inside(bracket_matcher_inside, e)
            if can_parse(v_1, tk_parser):
                # great, solved max information
                return v_1, 0
            else:
                # using inside not viable
                # replace bracket stuff and inside completey
                v_2 = remove_inside_group(bracket_matcher_inside, e)
                if can_parse(v_2, tk_parser):
                    # great, solved, just without bracket stuff
                    return v_2, 0
                else:
                    # might strip/subs required
                    v_3 = v_2
                    try:
                        v_3 = updated_texture_error_lookup[v_2]
                    except KeyError as ex:
                        pass
                    v_4 = strip_more(v_3)
                    if can_parse(v_4, t_parser):
                        # great, solved, just without bracket stuff
                        return v_4, 0
                    else:
                        logger.error(f'right inside branch: {t_e} {e} {v_2} {v_3} {v_4}')
                        return e, 1
                    
        else:
            # inside could NOT be parsed
            # replace bracket stuff and inside completey
            v_2 = remove_inside_group(bracket_matcher_inside, e)
            if can_parse(v_2, tk_parser):
                # great, solved, just without bracket stuff
                return v_2, 0
            else:
                # might strip/subs required
                v_3 = v_2
                try:
                    v_3 = updated_texture_error_lookup[v_2]
                except KeyError as ex:
                    pass
                v_4 = strip_more(v_3)
                if can_parse(v_4, tk_parser):
                    # great, solved, just without bracket stuff
                    return v_4, 0
                else:
                    logger.error(f'left inside branch: {t_e} {e} {v_2} {v_3} {v_4}')
                    return e, 1
    else:
        if can_parse(e, tk_parser):
            # great, solved, just without bracket stuff
            return e, 0
        else:
            # might strip/subs required
            v_4 = strip_more(e)
            if can_parse(v_4, tk_parser):
                # great, solved, just without bracket stuff
                return v_4, 0
            else:
                logger.error(f'outside branch: {t_e} {e} empty empty {v_4}')
                return e, 1


def consolidate_loimis_multiple_p(t_e: str, parsers_d: Dict) -> (str, int):
    e = t_e
    try:
        e = updated_texture_error_lookup[e]
    except KeyError as ex:
        pass
    v = bracket_matcher_inside.search(e)
    # test if has brackets
    if v is not None:
        t1 = v.group('inside')
        # test if inside can be parsed at all
        if can_parse_multiple(t1, parsers_d):
            # valid inside, test with inclusion in main string (just removing the brackets)
            v_1 = remove_brkt_inside(bracket_matcher_inside, e)
            if can_parse_multiple(v_1, parsers_d):
                # great, solved max information
                return v_1, 0
            else:
                # using inside not viable
                # replace bracket stuff and inside completey
                v_2 = remove_inside_group(bracket_matcher_inside, e)
                if can_parse_multiple(v_2, parsers_d):
                    # great, solved, just without bracket stuff
                    return v_2, 0
                else:
                    # might strip/subs required
                    v_3 = v_2
                    try:
                        v_3 = updated_texture_error_lookup[v_2]
                    except KeyError as ex:
                        pass
                    v_4 = strip_more(v_3)
                    if can_parse_multiple(v_4, parsers_d):
                        # great, solved, just without bracket stuff
                        return v_4, 0
                    else:
                        logger.error(f'right inside branch: {t_e} {e} {v_2} {v_3} {v_4}')
                        return e, 1
                    
        else:
            # inside could NOT be parsed
            # replace bracket stuff and inside completey
            v_2 = remove_inside_group(bracket_matcher_inside, e)
            if can_parse_multiple(v_2, parsers_d):
                # great, solved, just without bracket stuff
                return v_2, 0
            else:
                # might strip/subs required
                v_3 = v_2
                try:
                    v_3 = updated_texture_error_lookup[v_2]
                except KeyError as ex:
                    pass
                v_4 = strip_more(v_3)
                if can_parse_multiple(v_4, parsers_d):
                    # great, solved, just without bracket stuff
                    return v_4, 0
                else:
                    logger.error(f'left inside branch: {t_e} {e} {v_2} {v_3} {v_4}')
                    return e, 1
    else:
        if can_parse_multiple(e, parsers_d):
            # great, solved, just without bracket stuff
            return e, 0
        else:
            # might strip/subs required
            v_4 = strip_more(e)
            if can_parse_multiple(v_4, parsers_d):
                # great, solved, just without bracket stuff
                return v_4, 0
            else:
                logger.error(f'outside branch: {t_e} {e} empty empty {v_4}')
                return e, 1


def recursive_test_parse(word: str, pos: int, tracker: str, t_parser: Any, tk_parser: Any) -> str:
    if pos <= len(word):
        test_str = word[0:pos]
        if can_parse(test_str, t_parser) or can_parse(test_str, tk_parser):
            return recursive_test_parse(word, pos+1, test_str, t_parser, tk_parser)
        else:
            if pos < len(word):
                return recursive_test_parse(word, pos+1, tracker, t_parser, tk_parser)
            else:
                return tracker
    else:
        return tracker


def parse_test(l_str: List[str], t_parser: Any, tk_parser: Any) -> pd.Series:
    fixed_list = []
    has_errors = 0
    for e in l_str:
        try:
            v, err = consolidate_loimis(e, t_parser, tk_parser)
            fixed_list.append(v)
            has_errors += err
        except Exception as ex:
            logger.error(ex)
            logger.error(e)
            fixed_list.append('no_info')
            has_errors += err
    return pd.Series([fixed_list, has_errors])


def parse_test_dask(t_str: str, t_parser: Any, tk_parser: Any) -> pd.Series:
    err = 0
    fixed_list = []
    has_errors = 0
    for e in t_str.split('||'):
        try:
            v, err = consolidate_loimis(e, t_parser, tk_parser)
            fixed_list.append(v)
            has_errors += err
        except Exception as ex:
            logger.error(ex)
            logger.error(e)
            fixed_list.append('no_info')
            has_errors += err
    fixed_str = "||".join(fixed_list)
    return pd.Series((fixed_str, has_errors))


def parse_test_dask_multiple(t_str: str, parsers_d: Dict) -> pd.Series:
    err = 0
    fixed_list = []
    has_errors = 0
    for e in t_str.split('||'):
        try:
            v, err = consolidate_loimis_multiple_p(e, parsers_d)
            fixed_list.append(v)
            has_errors += err
        except Exception as ex:
            logger.error(ex)
            logger.error(e)
            fixed_list.append('no_info')
            has_errors += err
    fixed_str = "||".join(fixed_list)
    return pd.Series((fixed_str, has_errors))


def parse_reconstituate(l_str: List[str]) -> pd.Series:
    num_elems = len(l_str)
    has_no_info = 0
    result_loimis = ''
    if num_elems <= 0:
        return pd.Series(['no_info', 8])
    if num_elems > 1:
        for e in l_str:
            if 'no_info' in e:
                has_no_info += 1
        result_loimis = '/'.join(l_str)
        return pd.Series([result_loimis, has_no_info])
    if num_elems == 1:
        if 'no_info' in l_str[0]:
            has_no_info += 1
        result_loimis = l_str[0]
        return pd.Series([result_loimis, has_no_info])


def parse_reconstituate_dask(t_str: str) -> pd.Series:
    l_str = list(t_str.split("||"))
    num_elems = len(l_str)
    has_no_info = 0
    result_loimis = ''
    if num_elems <= 0:
        return pd.Series(('no_info', 8))
    if num_elems > 1:
        for e in l_str:
            if 'no_info' in e:
                has_no_info += 1
        result_loimis = '/'.join(l_str)
        return pd.Series((result_loimis, has_no_info))
    if num_elems == 1:
        if 'no_info' in l_str[0]:
            has_no_info += 1
        result_loimis = l_str[0]
        return pd.Series((result_loimis, has_no_info))


def load_default_texture_defensively(row: Dict[str, Any]) -> pd.Series:

    if row['loimis_reconst'] == 'no_info' or row['Loimis1'] is None:

        try:
            found = fillers_by_numbers.get(row['upd_siffer']).get('main_texture')

            return pd.Series([found, 0])
        except KeyError as ex:
            logger.error(ex)
            return pd.Series(['no_info', 1])
        except IndexError as ex:
            logger.error(ex)
            return pd.Series(['no_info', 1])
    else:
        return pd.Series([row['loimis_reconst'], 0])


def load_default_texture_defensively_dask(row: Dict[str, Any], fillers_by_numbers: Dict) -> pd.Series:

    if row['loimis_reconst'] == 'no_info' or row['Loimis1'] is None:

        try:
            found = fillers_by_numbers.get(row['upd_siffer']).get('main_texture')

            return pd.Series((found, 0))
        except KeyError as ex:
            logger.error(ex)
            return pd.Series(('no_info', 1))
        except IndexError as ex:
            logger.error(ex)
            return pd.Series(('no_info', 1))
        except AttributeError as ex:
            logger.error(ex)
            return pd.Series(('no_info', 1))
    else:
        return pd.Series((row['loimis_reconst'], 0))

##################################################################3
# HUMUS stuff
######################

def split_and_reduce(x):
    y = str(x)
    first_osa = y.replace(" ", ";").split(";")[0]
    return first_osa


def consolidate_num_bracket_humus(e: str) -> str:
    # r1 90(40) & r3 30-50(100)
    # r2 (70-100) & r4 50-90(40-50)
    # r2 (70-100) & r5 90(40-50)
    # r3 50-90(40-50) & r5 90(40-50)
    
    # r3 before r5
    r3 = bracket_matcher_num3.search(e)
    if r3 is not None:
        return repl_bt_reg(bracket_matcher_num3, e)
    
    # r4 before r2
    r4 = bracket_matcher_num4.search(e)
    if r4 is not None:
        return repl_bt_reg(bracket_matcher_num4, e)
        
    # r7 before r2
    r5 = bracket_matcher_num7.search(e)
    if r5 is not None:
        return repl_bt_reg(bracket_matcher_num7, e)
    
    # r3 before r1
    r1 = bracket_matcher_num1.search(e)
    if r1 is not None:
        return repl_bt_reg(bracket_matcher_num1, e)
    
    r2 = bracket_matcher_num2.search(e)
    if r2 is not None:
        return repl_bt_reg(bracket_matcher_num2, e)
    
    r6 = bracket_matcher_num6.search(e)
    if r6 is not None:
        return repl_bt_reg(bracket_matcher_num6, e)
    
    return e


def test_brackets_huumus_dask(l_str: str) -> str:
    fixed_list = []
    for e in l_str.split('/'):
        txt = subscripts_consolidate(e)
        result = bracket_matcher.search(txt)
        if result is not None:
            num_fix = consolidate_num_bracket_humus(txt)
            # if num_fix == txt:
            #     global_brackets_list.append(txt)
            fixed_list.append(num_fix)
        else:
            fixed_list.append(txt)
    return "/".join(fixed_list)


def cleanup_replace(x: str) -> str:
    a1 = x.replace("+","-")
    return a1


def strip_more_humus(x: str) -> str:
    no_comma = x.split(",")
    bm = open_bracket_matcher.search(no_comma[0])
    nxt = ''
    if bm is not None:
        nxt = remove_inside_group(open_bracket_matcher, no_comma[0])
    else:
        nxt = no_comma[0]
    xdash = cleanup_replace(nxt)
    return xdash


def huumus_grammar():

    def amp1(): return '₁'

    def amp2(): return '₂'

    def amp3(): return '₃'

    def amp4(): return '₄'

    def amp5(): return '₅'

    def amplifiers(): return [amp1, amp2, amp3, amp4, amp5]


    def depth_number(): return RegExMatch(r'\d+')  # at least two digits for depth numbers (well might as well be one)

    def kodu_depth_number(): return depth_number, Optional(amplifiers)

    def depth_range(): return depth_number, ZeroOrMore('-', depth_number)

    def th(): return 'th'  # th 15 või th 15-20 toorhuumusliku horisondi tüsedus, en: raw humus thickness

    # def t(): return Combine( [ (Not(tls), 't'), (Not(tsl), 't')] )
    def t(): return 't'

    def th_turfs(): return th, depth_range

    def t_turfs(): return t, Optional(amplifiers), depth_range

    def huumus_hor(): return depth_range

    def kodu_depth_range(): return kodu_depth_number, ZeroOrMore('-', kodu_depth_number)

    
    def th_turfs_grammar(): return th_turfs, EOF

    def t_turfs_grammar(): return t_turfs, EOF

    def huumus_hor_grammar(): return huumus_hor, EOF
    
    def kodu_depth_range_grammar(): return kodu_depth_range, EOF
    
    th_turfs_parser = ParserPython(th_turfs_grammar,memoization=False)
    t_turfs_parser = ParserPython(t_turfs_grammar,memoization=False)
    huumus_hor_parser = ParserPython(huumus_hor_grammar,memoization=False)
    kodu_depth_range_parser = ParserPython(kodu_depth_range_grammar,memoization=False)

    parsers_d = {'th_turfs_parser': th_turfs_parser,
               't_turfs_parser':t_turfs_parser,
               'huumus_hor_parser':huumus_hor_parser,
               'kodu_depth_range_parser':kodu_depth_range_parser }

    return parsers_d


def initial_huumus_clean(t_e: str, parsers_d: Dict) -> (str, int):
    e = t_e
    try:
        e = updated_texture_error_lookup[e]
    except KeyError as ex:
        pass
    v = bracket_matcher_inside.search(e)
    # test if has brackets
    if v is not None:
        t1 = v.group('inside')
        # test if inside can be parsed at all
        if can_parse_multiple(t1, parsers_d):
            # valid inside, test with inclusion in main string (just removing the brackets)
            v_1 = remove_brkt_inside(bracket_matcher_inside, e)
            if can_parse_multiple(v_1, parsers_d):
                # great, solved max information
                return v_1, 0
            else:
                # using inside not viable
                # replace bracket stuff and inside completey
                v_2 = remove_inside_group(bracket_matcher_inside, e)
                if can_parse_multiple(v_2, parsers_d):
                    # great, solved, just without bracket stuff
                    return v_2, 0
                else:
                    # might strip/subs required
                    v_3 = v_2
                    try:
                        v_3 = updated_texture_error_lookup[v_2]
                    except KeyError as ex:
                        pass
                    v_4 = strip_more_humus(v_3)
                    if can_parse_multiple(v_4, parsers_d):
                        # great, solved, just without bracket stuff
                        return v_4, 0
                    else:
                        logger.error(f'right inside branch: {t_e} {e} {v_2} {v_3} {v_4}')
                        return e, 1
                    
        else:
            # inside could NOT be parsed
            # replace bracket stuff and inside completey
            v_2 = remove_inside_group(bracket_matcher_inside, e)
            if can_parse_multiple(v_2, parsers_d):
                # great, solved, just without bracket stuff
                return v_2, 0
            else:
                # might strip/subs required
                v_3 = v_2
                try:
                    v_3 = updated_texture_error_lookup[v_2]
                except KeyError as ex:
                    pass
                v_4 = strip_more_humus(v_3)
                if can_parse_multiple(v_4, parsers_d):
                    # great, solved, just without bracket stuff
                    return v_4, 0
                else:
                    logger.error(f'left inside branch: {t_e} {e} {v_2} {v_3} {v_4}')
                    return e, 1
    else:
        if can_parse_multiple(e, parsers_d):
            # great, solved, just without bracket stuff
            return e, 0
        else:
            # might strip/subs required
            v_4 = strip_more_humus(e)
            if can_parse_multiple(v_4, parsers_d):
                # great, solved, just without bracket stuff
                return v_4, 0
            else:
                logger.error(f'outside branch: {t_e} {e} empty empty {v_4}')
                return e, 1


def parse_test_humus_multiple(t_str: str, parsers_d: Dict) -> pd.Series:
    err = 0
    fixed_list = []
    has_errors = 0
    for e in t_str.split('/'):
        try:
            v, err = initial_huumus_clean(e, parsers_d)
            fixed_list.append(v)
            has_errors += err
        except Exception as ex:
            logger.error(ex)
            logger.error(e)
            fixed_list.append('no_info')
            has_errors += err
    fixed_str = "/".join(fixed_list)
    return pd.Series((fixed_str, has_errors))

