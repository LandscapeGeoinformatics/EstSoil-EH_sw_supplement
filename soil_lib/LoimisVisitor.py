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

import numpy as np
import pandas as pd

import traceback
import logging

from arpeggio import PTNodeVisitor, visit_parse_tree, ParserPython

from .LoimisLookups import updated_texture_error_lookup, texture_rules, kores_amp_strength

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


##############################
# Arpeggio-based Grammar for the Texture value encoding (loimis)
#################################

###############################################################################################
# NEW NEW NEW
################################################################################

class LpVisitorV2(PTNodeVisitor):

    def visit_loimisp_short_summary(self, node, children):
        """
        # new
        def loimisp_short(): return constituent, EOF
        """
        if self.debug:
            logger.debug("loimisp_short node {}".format(node.value))
            logger.debug("loimisp_short({}) {}".format(len(children), children))

            if len(children) > 2:
                logger.info("loimisp_short untypically long({}) {}".format(
                    len(children), children))

        tmp = [p for p in children if isinstance(p, dict)]

        # we just reduce to one
        loimis_def = {"type": "loimis",
                      "count": len(tmp),
                      "soilparts": tmp}
        
        return loimis_def


    # TODO
    def visit_loimisp_short_dk_1(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_2(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_3(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_4(self, node, children):
        return self.visit_loimisp_short_summary(node, children)

    def visit_loimisp_short_dk_1a(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_2a(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_3a(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_4a(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_5(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    def visit_loimisp_short_dk_6(self, node, children):
        return self.visit_loimisp_short_summary(node, children)
    
    
    def visit_constituent(self, node, children):
        if self.debug:
            logger.debug("constituent node {}".format(node.value))
            logger.debug("constituent({}) {}".format(len(children), children))

            if len(children) > 2:
                logger.info("constituent untypically long({}) {}".format(
                    len(children), children))

        # one compound element i.e. one type combo of kores, turf,soil with amps or karbonated and depth,
        if len(children) > 2:
            logger.info("constituent untypically long ({}) {}".format(
                len(children), children))
        # we reduce to one skeleton (koores) and one peat (turvas)/fine earth (peenes)?
        # typically kores-peenes, turf-peenes, sometime kores-kores-(turf/peenes), then cut one kores off
        # basically check for repetations (+1 same?) and pop the following, keep first

        count_kores = len([k for k in children if k['type'] == 'kores'])
        count_turfs = len([k for k in children if k['type'] == 'turfs'])
        count_peenes = len([k for k in children if k['type'] == 'peenes'])

        has_kores = False
        has_turfs = False
        has_peenes = False

        constituents_final = []

        for x in children:
            if isinstance(x, dict) and x['type'] == 'kores' and has_kores == False:
                has_kores = True
                constituents_final.append(x)
            elif isinstance(x, dict) and x['type'] == 'kores' and has_kores == True:
                logger.debug('dropping a kores, because of double')
        
        for x in children:
            if isinstance(x, dict) and x['type'] == 'turfs' and has_peenes == False:
                has_peenes = True
                constituents_final.append(x)
            elif isinstance(x, dict) and x['type'] == 'turfs' and has_peenes == True:
                logger.debug('dropping a turfs, because of double')
        
        for x in children:
            if isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == False:
                has_peenes = True
                constituents_final.append(x)
            elif isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == True:
                logger.debug('dropping a peenes, because of double')

        return {"count": len(constituents_final),
                "constituents": constituents_final}
    

    # TODO
    def visit_double_kores_consitutent(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_1(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_2(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_3(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_4(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_1a(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_2a(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_3a(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_4a(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_5(self, node, children):
        return self.visit_constituent(node, children)
    
    def visit_double_kores_consitutent_6(self, node, children):
        return self.visit_constituent(node, children)


    def visit_kores(self, node, children):
        if self.debug:
            logger.debug("kores node {}".format(node.value))
            logger.debug("kores ({}) {}".format(len(children), children))

        # is sum list of kores variants, with amp, k-plus etc
        dat = {"type": "kores"}
        t_child = children

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)
        else:
            dat.update({"karbonaat": False})
        
        
        if len(t_child) > 0 and isinstance(t_child[0], str):
            dat.update({"code": t_child[0]})
            t_child.pop(0)
            
        # an optional amp for the now scattered/hierarchical kores
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'amp' in t_child[0]:
            dat.update(t_child[0])
            t_child.pop(0)
        else:
            dat.update({"amp": False})
        
        try:
            logger.info(f"at that place: {t_child} {type(t_child)}  {type(t_child[0])}  {type(t_child[0][1])}")
            sub_amp = t_child[0][1]
            if isinstance(sub_amp, dict) and 'amp' in sub_amp:
                dat.update(sub_amp)
                dat.update({"code": t_child[0][0]})
                t_child.pop(0)
        except KeyError:
            pass
        except IndexError:
            pass

        # an optional depth related
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'depth' in t_child[0]:
            dat.update(t_child[0])
        else:
            dat.update({"depth": False})

        return dat
    
    
    def visit_amp_kores(self, node, children):
        return self.visit_kores(node, children)
    
    def visit_no_amp_kores(self, node, children):
        return self.visit_kores(node, children)
    
    def visit_p_kores(self, node, children):
        return self.visit_kores(node, children)
    
    def visit_p(self, node, children):
        return node.value
    
    def visit_skeleton_for_amp(self, node, children):
        if self.debug:
            logger.debug("skeleton_with_amp node {}".format(node.value))
            logger.debug("skeleton_with_amp {}".format(children))
        # contains single kores variants from list
        if len(children) == 1:
            return children[0]
        if len(children) == 2:
            logger.info(f"children kores with amp: {children}")
            return children
        else:
            logger.info("warning too many skeleton_with_amp here {}".format(children))
            return children


    def visit_skeleton_no_amp(self, node, children):
        if self.debug:
            logger.debug("skeleton_no_amp node {}".format(node.value))
            logger.debug("skeleton_no_amp {}".format(children))
        # contains single kores variants from list
        if len(children) == 1:
            return children[0]
        else:
            logger.info("warning too many visit_skeleton_no_amp here {}".format(children))
            return children


    def visit_peenes(self, node, children):
        if self.debug:
            logger.debug("peenes node {}".format(node.value))
            logger.debug("peenes ({}) {}".format(len(children), children))

        # is sum list of peenes variants incl amps, depth
        dat = {"type": "peenes"}
        t_child = children

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)
        else:
            dat.update({"karbonaat": False})

        # the compulsory peenes code
        dat.update({"code": children[0]})
        t_child.pop(0)

        # an optional amp for the peenes
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'amp' in t_child[0]:
            dat.update(t_child[0])
            t_child.pop(0)
        else:
            dat.update({"amp": False})

        # the kPlus could b on this position, too
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)

        # an optional depth related
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'depth' in t_child[0]:
            dat.update(t_child[0])
        else:
            dat.update({"depth": False})

        return dat

    
    def visit_peenes_t(self, node, children):
        return self.visit_peenes(node, children)
    
    def visit_peenes_p(self, node, children):
        return self.visit_peenes(node, children)
    
    
    def visit_peenes_list(self, node, children):
        """
        def peenes_list(): return [plsl, pl, tsl, tls, dk, sl, ls, s, l]
        """
        if self.debug:
            logger.debug("peenes_list node {}".format(node.value))
            logger.debug("peenes_list {}".format(children))
        # contains single peenes variants from list
        if len(children) == 1:
            return children[0]
        else:
            logger.info(
                "warning too many peenes_list here {}".format(children))
            return children
    
    
    def visit_peenes_list_t(self, node, children):
        return self.visit_peenes_list(node, children)
    
    def visit_peenes_list_p(self, node, children):
        return self.visit_peenes_list(node, children)


    def visit_turfs(self, node, children):
        if self.debug:
            logger.debug("turfs node {}".format(node.value))
            logger.debug("turfs {}".format(children))

        # contains turfs variants incl amps, depth
        dat = {"type": "turfs"}
        t_child = children
        
        print('node')
        print(node)
        print('children')
        print(children)
        dat.update({"code": node.value})

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)
        else:
            dat.update({"karbonaat": False})

        # the compulsory peenes/turvas code
        dat.update({"code": children[0]})
        t_child.pop(0)

        # an optional amp for the peenes
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'amp' in t_child[0]:
            dat.update(t_child[0])
            t_child.pop(0)
        else:
            dat.update({"amp": False})

        # an optional depth related
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'depth' in t_child[0]:
            dat.update(t_child[0])
        else:
            dat.update({"depth": False})

        return dat

    
    def visit_th_turfs(self, node, children):
        return self.visit_turfs(node, children)
    
    def visit_t_turfs(self, node, children):
        return self.visit_turfs(node, children)
    
    
    def visit_t(self, node, children):
        return node.value
    
    def visit_th(self, node, children):
        return node.value
        
    def visit_amplifiers(self, node, children):
        """
        def amplifiers(): return [amp1, amp2, amp3, amp4, amp5]
        """
        if self.debug:
            logger.debug("amplifiers node {}".format(node.value))
            logger.debug("amplifiers {}".format(children))

        # compound list of amps, should return only one number 1-5, might not be needed for texture or similar?
        if node.value == '₁':
            return {"amp": 1}
        if node.value == '₂':
            return {"amp": 2}
        if node.value == '₃':
            return {"amp": 3}
        if node.value == '₄':
            return {"amp": 4}
        if node.value == '₅':
            return {"amp": 5}
        return {"amp": int(node.value)}


    def visit_depth_range(self, node, children):
        """
        def depth_range(): return Optional(kPlus), depth_number, ZeroOrMore('-', depth_number)
        """
        if self.debug:
            logger.debug("depth_range {}".format(children))
            logger.debug("depth_range node {}".format(node.value))

        dat = {"depth": {}}
        t_child = children

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat['depth'].update({"deeper_than": True})
            t_child.pop(0)

        if len(t_child) == 2:
            logger.info("depth is from to {}".format(children))
            dat['depth'].update(
                {"range": True,
                 "from": t_child[0],
                 "to": t_child[1]
                 })

        if len(t_child) == 1:
            dat['depth'].update(
                {"range": False,
                 "to": t_child[0]
                 })
        # compound operator for depth elements, either has one number, or two numbers separated with '-'
        # should return the calculated depth horizon?
        return dat


    def visit_depth_number(self, node, children):
        """
        def depth_number(): return RegExMatch(r'\\d\\d+') # at least two digits for depth numbers
        """
        if self.debug:
            logger.debug("depth_number {}".format(children))
        # itself might not be needed if handled already by range?
        # here we do the SWAT cm to mm conversion
        return float(node.value)*10


    def visit_kPlus(self, node, children):
        if self.debug:
            logger.debug("kPlus {}".format(children))
        # karbonated option, should either return yes or no? might not be needed?
        return node.value


class LpVisitor(PTNodeVisitor):
    """The main Arpeggio Visitor implementaion to evaluate the found expressions and build a nice data structure
    
    Arguments:
        PTNodeVisitor -- the class to be inherited from Arpeggio

    Returns:
        Visitor
    """
    def visit_soilParts(self, node, children):
        """
        # new
        def soilParts(): return constituent, ZeroOrMore(vertiSep, constituent), EOF

        # old
        def soilParts(): return paraSeq, ZeroOrMore(vertiSep, paraSeq)
        """
        if self.debug:
            logger.debug("soilParts node {}".format(node.value))
            logger.debug("soilParts({}) {}".format(len(children), children))
        # soil parts is the sum list of the sequence of constituents,
        # separated by a vertiSep (at different depths)

        tmp = [p for p in children if isinstance(p, dict)]

        # we just reduce to one
        loimis_def = {"type": "loimis",
                      "count": len(tmp),
                      "soilparts": tmp}
        
        return loimis_def
    

    # TODO
    def visit_soilParts_dk(self, node, children):
        """
        # new
        soilParts_dk(): return double_kores_consitutent, ZeroOrMore(vertiSep, double_kores_consitutent), EOF
        """
        return self.visit_soilParts(node, children)
    

    def visit_loimisp_short(self, node, children):
        """
        # new
        def loimisp_short(): return constituent, EOF
        """
        if self.debug:
            logger.debug("loimisp_short node {}".format(node.value))
            logger.debug("loimisp_short({}) {}".format(len(children), children))

            if len(children) > 2:
                logger.info("loimisp_short untypically long({}) {}".format(
                    len(children), children))

        tmp = [p for p in children if isinstance(p, dict)]
        return tmp[0]


    # TODO
    def visit_loimisp_short_dk(self, node, children):
        """
        # new
        def loimisp_short_dk(): return double_kores_consitutent, EOF
        """
        return self.visit_loimisp_short(node, children)
    

    def visit_constituent(self, node, children):
        """
        # new
        def constituent(): return Optional(kores), Optional(turfs), Optional(peenes)

        # old
        def constituent(): return ZeroOrMore(kores), ZeroOrMore([peenes, turfs, Optional(alternateComma)])
        """
        if self.debug:
            logger.debug("constituent node {}".format(node.value))
            logger.debug("constituent({}) {}".format(len(children), children))

            if len(children) > 2:
                logger.info("constituent untypically long({}) {}".format(
                    len(children), children))

        # one compound element i.e. one type combo of kores, turf,soil with amps or karbonated and depth,
        if len(children) > 2:
            logger.info("constituent untypically long ({}) {}".format(
                len(children), children))
        # we reduce to one skeleton (koores) and one peat (turvas)/fine earth (peenes)?
        # typically kores-peenes, turf-peenes, sometime kores-kores-(turf/peenes), then cut one kores off
        # basically check for repetations (+1 same?) and pop the following, keep first

        count_kores = len([k for k in children if k['type'] == 'kores'])
        count_turfs = len([k for k in children if k['type'] == 'turfs'])
        count_peenes = len([k for k in children if k['type'] == 'peenes'])

        has_kores = False
        has_turfs = False
        has_peenes = False

        constituents_final = []

        for x in children:
            if isinstance(x, dict) and x['type'] == 'kores' and has_kores == False:
                has_kores = True
                constituents_final.append(x)
            elif isinstance(x, dict) and x['type'] == 'kores' and has_kores == True:
                logger.debug('dropping a kores, because of double')
        
        for x in children:
            if isinstance(x, dict) and x['type'] == 'turfs' and has_peenes == False:
                has_peenes = True
                constituents_final.append(x)
            elif isinstance(x, dict) and x['type'] == 'turfs' and has_peenes == True:
                logger.debug('dropping a turfs, because of double')
        
        for x in children:
            if isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == False:
                has_peenes = True
                constituents_final.append(x)
            elif isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == True:
                logger.debug('dropping a peenes, because of double')

        return {"count": len(constituents_final),
                "constituents": constituents_final}
    

    # TODO
    def visit_double_kores_consitutent(self, node, children):
        """
        # new
        def double_kores_consitutent(): return ZeroOrMore(kores), Optional(turfs), Optional(peenes)
        """
        return self.visit_constituent(node, children)


    def visit_kores(self, node, children):
        """
        # new
        def kores(): return Optional(kPlus), [ skeleton_no_amp, skeleton_with_amp ], Optional(depth_range)
        
        # old
        def kores(): return Optional(kPlus), kores_list, Optional(amplifiers), Optional(depth_range)
        """
        if self.debug:
            logger.debug("kores node {}".format(node.value))
            logger.debug("kores ({}) {}".format(len(children), children))

        # is sum list of kores variants, with amp, k-plus etc
        dat = {"type": "kores"}
        t_child = children

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)
        else:
            dat.update({"karbonaat": False})
        
        
        if len(t_child) > 0 and isinstance(t_child[0], str):
            dat.update({"code": t_child[0]})
            t_child.pop(0)
            
        # an optional amp for the now scattered/hierarchical kores
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'amp' in t_child[0]:
            dat.update(t_child[0])
            t_child.pop(0)
        else:
            dat.update({"amp": False})
        
        try:
            logger.info(f"at that place: {t_child} {type(t_child)}  {type(t_child[0])}  {type(t_child[0][1])}")
            sub_amp = t_child[0][1]
            if isinstance(sub_amp, dict) and 'amp' in sub_amp:
                dat.update(sub_amp)
                dat.update({"code": t_child[0][0]})
                t_child.pop(0)
        except KeyError:
            pass
        except IndexError:
            pass

        # an optional depth related
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'depth' in t_child[0]:
            dat.update(t_child[0])
        else:
            dat.update({"depth": False})

        return dat
    

    def visit_skeleton_with_amp(self, node, children):
        """
        def skeleton_with_amp(): return [  r_norm, v_norm, kb_norm, k_norm ], Optional(amplifiers)
        """
        if self.debug:
            logger.debug("skeleton_with_amp node {}".format(node.value))
            logger.debug("skeleton_with_amp {}".format(children))
        # contains single kores variants from list
        if len(children) == 1:
            return children[0]
        if len(children) == 2:
            logger.info(f"children kores with amp: {children}")
            return children
        else:
            logger.info("warning too many skeleton_with_amp here {}".format(children))
            return children


    def visit_skeleton_no_amp(self, node, children):
        """
        def skeleton_no_amp(): return [  no_info, pk, kr, p, d, lu, ck ]
        """
        if self.debug:
            logger.debug("skeleton_no_amp node {}".format(node.value))
            logger.debug("skeleton_no_amp {}".format(children))
        # contains single kores variants from list
        if len(children) == 1:
            return children[0]
        else:
            logger.info("warning too many visit_skeleton_no_amp here {}".format(children))
            return children


    def visit_peenes(self, node, children):
        """
        def peenes(): return Optional(kPlus), peenes_list, Optional(amplifiers), Optional(kPlus), Optional(depth_range)
        """
        if self.debug:
            logger.debug("peenes node {}".format(node.value))
            logger.debug("peenes ({}) {}".format(len(children), children))

        # is sum list of peenes variants incl amps, depth
        dat = {"type": "peenes"}
        t_child = children

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)
        else:
            dat.update({"karbonaat": False})

        # the compulsory peenes code
        dat.update({"code": children[0]})
        t_child.pop(0)

        # an optional amp for the peenes
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'amp' in t_child[0]:
            dat.update(t_child[0])
            t_child.pop(0)
        else:
            dat.update({"amp": False})

        # the kPlus could b on this position, too
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)

        # an optional depth related
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'depth' in t_child[0]:
            dat.update(t_child[0])
        else:
            dat.update({"depth": False})

        return dat


    def visit_peenes_list(self, node, children):
        """
        def peenes_list(): return [plsl, pl, tsl, tls, dk, sl, ls, s, l]
        """
        if self.debug:
            logger.debug("peenes_list node {}".format(node.value))
            logger.debug("peenes_list {}".format(children))
        # contains single peenes variants from list
        if len(children) == 1:
            return children[0]
        else:
            logger.info(
                "warning too many peenes_list here {}".format(children))
            return children


    def visit_turfs(self, node, children):
        """
        def turfs(): return Optional(kPlus), [th, t], Optional(amplifiers), Optional(depth_range)
        """
        if self.debug:
            logger.debug("turfs node {}".format(node.value))
            logger.debug("turfs {}".format(children))

        # contains turfs variants incl amps, depth
        dat = {"type": "turfs"}
        t_child = children

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat.update({"karbonaat": True})
            t_child.pop(0)
        else:
            dat.update({"karbonaat": False})

        # the compulsory peenes/turvas code
        dat.update({"code": children[0]})
        t_child.pop(0)

        # an optional amp for the peenes
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'amp' in t_child[0]:
            dat.update(t_child[0])
            t_child.pop(0)
        else:
            dat.update({"amp": False})

        # an optional depth related
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'depth' in t_child[0]:
            dat.update(t_child[0])
        else:
            dat.update({"depth": False})

        return dat


    def visit_amplifiers(self, node, children):
        """
        def amplifiers(): return [amp1, amp2, amp3, amp4, amp5]
        """
        if self.debug:
            logger.debug("amplifiers node {}".format(node.value))
            logger.debug("amplifiers {}".format(children))

        # compound list of amps, should return only one number 1-5, might not be needed for texture or similar?
        if node.value == '₁':
            return {"amp": 1}
        if node.value == '₂':
            return {"amp": 2}
        if node.value == '₃':
            return {"amp": 3}
        if node.value == '₄':
            return {"amp": 4}
        if node.value == '₅':
            return {"amp": 5}
        return {"amp": int(node.value)}


    def visit_depth_range(self, node, children):
        """
        def depth_range(): return Optional(kPlus), depth_number, ZeroOrMore('-', depth_number)
        """
        if self.debug:
            logger.debug("depth_range {}".format(children))
            logger.debug("depth_range node {}".format(node.value))

        dat = {"depth": {}}
        t_child = children

        # optional first kPlus
        if len(t_child) > 0 and t_child[0] == "+":
            dat['depth'].update({"deeper_than": True})
            t_child.pop(0)

        if len(t_child) == 2:
            logger.info("depth is from to {}".format(children))
            dat['depth'].update(
                {"range": True,
                 "from": t_child[0],
                 "to": t_child[1]
                 })

        if len(t_child) == 1:
            dat['depth'].update(
                {"range": False,
                 "to": t_child[0]
                 })
        # compound operator for depth elements, either has one number, or two numbers separated with '-'
        # should return the calculated depth horizon?
        return dat


    def visit_depth_number(self, node, children):
        """
        def depth_number(): return RegExMatch(r'\\d\\d+') # at least two digits for depth numbers
        """
        if self.debug:
            logger.debug("depth_number {}".format(children))
        # itself might not be needed if handled already by range?
        # here we do the SWAT cm to mm conversion
        return float(node.value)*10


    def visit_kPlus(self, node, children):
        if self.debug:
            logger.debug("kPlus {}".format(children))
        # karbonated option, should either return yes or no? might not be needed?
        return node.value


    def visit_vertiSep(self, node, children):
        if self.debug:
            logger.debug("vertiSep {}".format(children))
        # vertical separator, should be known at appropriate higher place for depth calculations
        return node.value


######################################################
#
#####################################################

def can_parse_multiple_get_parser(e: str, parsers_d: Dict) -> Any:
    for parser_name in parsers_d.keys():
        parser = parsers_d[parser_name]
        done = False
        try:
            parse_tree = parser.parse(e)
            if parse_tree is not None:
                result_dict = visit_parse_tree(parse_tree, LpVisitor(debug=False))
                done = True
                return parser
        except Exception as ex:
            continue

        if done:
            break


def loimis_grammar_product(input_expr: str, parser: ParserPython) -> pd.Series:
    """applies the grammar to a textue string and builds a usable data structure

    Arguments:
        input_expr {str} -- the texture string that will be evaluated by the grammar
        parser {ParserPython} -- the Parser object that has the grammar definitions loaded

    Returns:
        pd.Series -- a Series with the defined field for adding to the soil map DataFrame
    """
    return_empty = {'type': 'loimis',
                    'count': 0,
                    'message': 'empty_loimis',
                    'soilparts': []
                    }
    return_error = {'type': 'loimis',
                    'count': 0,
                    'message': 'parse_error',
                    'soilparts': []
                    }
    try:
        if input_expr == None or input_expr == 'no_info':
            parse_info = return_empty['message']
            return_empty['code'] = input_expr
            return pd.Series([return_empty, parse_info])
        parse_tree = parser.parse(input_expr)
        # logger.debug("expr: {} >> parse_tree: {}".format(input_expr, parse_tree))
        result_dict = visit_parse_tree(parse_tree, LpVisitor(debug=False))
        result_dict['code'] = input_expr
        # logger.info("expr: {} >> result_dict: {}".format(input_expr, result_dict))
        result_dict['message'] = 'successful'
        return pd.Series([result_dict, input_expr])
    except Exception as ex:
        logger.error(ex)
        logger.error("error expr >> {}".format(input_expr))
        return_error['code'] = input_expr
        parse_info = return_error['message']
        return pd.Series([return_error, parse_info])


def loimis_grammar_product_dask(input_expr: str, parser: ParserPython) -> pd.Series:
    """applies the grammar to a textue string and builds a usable data structure

    Arguments:
        input_expr {str} -- the texture string that will be evaluated by the grammar
        parser {ParserPython} -- the Parser object that has the grammar definitions loaded

    Returns:
        pd.Series -- a Series with the defined field for adding to the soil map DataFrame
    """
    return_empty = {'type': 'loimis',
                    'count': 0,
                    'message': 'empty_loimis',
                    'soilparts': []
                    }
    return_error = {'type': 'loimis',
                    'count': 0,
                    'message': 'parse_error',
                    'soilparts': []
                    }
    try:
        if input_expr == None or input_expr == 'no_info':
            parse_info = return_empty['message']
            return_empty['code'] = input_expr
            return pd.Series((return_empty, parse_info))
        parse_tree = parser.parse(input_expr)
        # logger.debug("expr: {} >> parse_tree: {}".format(input_expr, parse_tree))
        result_dict = visit_parse_tree(parse_tree, LpVisitor(debug=False))
        result_dict['code'] = input_expr
        # logger.info("expr: {} >> result_dict: {}".format(input_expr, result_dict))
        result_dict['message'] = 'successful'
        return pd.Series((result_dict, input_expr))
    except Exception as ex:
        logger.error(ex)
        logger.error("error expr >> {}".format(input_expr))
        return_error['code'] = input_expr
        parse_info = return_error['message']
        return pd.Series((return_error, parse_info))


def loimis_grammar_product_dask_multiple(input_expr: str, parsers_d: Dict) -> pd.Series:
    """applies the grammar to a textue string and builds a usable data structure"""
    
    return_base = {'type': 'loimis',
                    'count': 0,
                    'message': 'empty_loimis',
                    'soilparts': [],
                    'code': input_expr
                    }
    return_error = {'type': 'loimis',
                    'count': 0,
                    'message': 'parse_error',
                    'soilparts': [],
                    'code': input_expr
                    }
    if input_expr == None or input_expr == 'no_info':
        parse_info = return_base['message']
        return pd.Series((return_base, parse_info))
    
    return_soilparts_list = []
    
    for e in input_expr.split('/'):
        result_dict = {'count': 0, 'constituents': []}
        try:
            parser = can_parse_multiple_get_parser(e, parsers_d)
            if parser is not None:
                parse_tree = parser.parse(e)
                # logger.debug("expr: {} >> parse_tree: {}".format(e, parse_tree))
                result_parsed = visit_parse_tree(parse_tree, LpVisitorV2(debug=False))
                
                if result_parsed['count'] > 0:
                    result_dict['constituents_code'] = e
                    const_list = result_parsed['soilparts'][0]['constituents']
                    result_dict['count'] = len(const_list)
                    # logger.info("expr: {} >> result_dict: {}".format(e, result_dict))
                    result_dict['message'] = 'successful'
                    result_dict['constituents'] = const_list
                    return_soilparts_list.append(result_dict)
                    
                elif result_parsed['count'] > 1:
                    logger.error("error expr more than one part >> {}".format(result_parsed))
                else:
                    logger.error("error expr0 parts >> {}".format(result_parsed))
            else:
                logger.error("error expr >> {}".format(e))
                parse_info = return_error['message']
                return_error['message'] = "parser is None, cant be parsed"
                return pd.Series((return_error, parse_info))
        except Exception as ex:
            logger.error(ex)
            logger.error("error expr >> {}".format(e))
            parse_info = return_error['message']
            return_error['message'] = ex
            return pd.Series((return_error, parse_info))
    
    return_base['count'] = len(return_soilparts_list)
    return_base['message'] = 'successful'
    return_base['soilparts'] = return_soilparts_list
    return pd.Series((return_base, input_expr))


def test_layer_depths(grammar_loimis: Dict[str, Any]) -> pd.Series:
    """examines the grammar data structure and derives the layer information
    Arguments:
        grammar_loimis {Dict[str, Any]} -- the Pandas column with the soil map and the grammar data structure
    Returns:
        pd.Series -- a Series with the defined field for adding to the soil map DataFrame
    """
    start_depth = 0
    default_depth = 1000
    return_code = {'layers': 1, 'depths': [default_depth]}
    depth_builder = []
    nlayers = 0
    SOL_ZMX = 0
    SOL_Z1 = 0
    SOL_Z2 = 0
    SOL_Z3 = 0
    SOL_Z4 = 0
    try:
        if isinstance(grammar_loimis, dict) and grammar_loimis.get('count', 0) <= 0:
            # logger.info('problem?')
            # return return_code
            pd.Series([nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4])
        # logger.info(result)
        # logger.info(type(result))
        soilparts_count = grammar_loimis['count']
        return_code.update({'layers': int(soilparts_count)})
        for y in range(0, soilparts_count):
            constituents_obj = grammar_loimis['soilparts'][y]
            constituents_count = constituents_obj['count']
            constituents_arr = constituents_obj['constituents']
            # one soil part can have multiple constituents e.g. r and ls, but usually only one depth
            depth_for_soilpart: Dict = {}
            depth_found_for_soilpart = False
            range_val = False
            from_val = 0
            to_val = 0
            depth_avg = 0
            for x in range(0, constituents_count):
                constitu = constituents_arr[x]
                if constitu['depth']:
                    depth_for_soilpart = constitu['depth']
                    depth_found_for_soilpart = True
            if depth_found_for_soilpart:
                try:
                    range_val = depth_for_soilpart['range']
                    # if 'deeper_than' in depth_for_soilpart:
                    #     deeper_than_val = depth_for_soilpart['deeper_than']
                    if 'from' in depth_for_soilpart:
                        from_val = int(depth_for_soilpart['from'])
                    if 'to' in depth_for_soilpart:
                        to_val = int(depth_for_soilpart['to'])
                    if range_val:
                        depth_avg = int((to_val + from_val)/2)
                    else:
                        depth_avg = to_val
                except KeyError as ex:
                    tb = traceback.format_exc()
                    logger.info("keyerror at {}".format(ex))
                    logger.info(tb)
            else:
                down_to_def = default_depth - start_depth
                depth_avg = down_to_def
            thickness = start_depth+depth_avg - start_depth
            if thickness <= 0:
                # logger.info('layers thickness of layer/soilpart {} is zero( thickness{},depth_avg{},start_depth{} )'
                #       .format(y+1, thickness,depth_avg,start_depth))
                return_code.update({'layers': int(soilparts_count-1)})
            else:
                if y == 0:
                    return_code.update({'depths': []})
                depth_builder.append(int(thickness))
                start_depth = start_depth+depth_avg
        return_code.update({'depths': depth_builder})
        # result = return_code
        nlayers = return_code['layers']
        depths = return_code['depths']
        SOL_ZMX = sum(depths)
        SOL_Z1 = depths[0]
        if len(depths) > 1 and nlayers > 1:
            SOL_Z2 = depths[1]
        if len(depths) > 2 and nlayers > 2:
            SOL_Z3 = depths[2]
        if len(depths) > 3 and nlayers > 3:
            SOL_Z4 = depths[3]
        return pd.Series([nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4])
    except KeyError as ex:
        tb = traceback.format_exc()
        logger.info("keyerror at {}".format(ex))
        logger.info(tb)
        return pd.Series([nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4])
    except Exception as ex:
        logger.info(ex)
        return pd.Series([nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4])


def has_only_no_info_constituents(constituents_arr):
    only_no_info = True
    for elem in constituents_arr:
        if elem.get('code', None) is not None and not elem.get('code', 'should_never') == 'no_info':
            only_no_info = False
            return only_no_info
    return only_no_info

        
def test_layer_depths_dask(grammar_loimis: Dict[str, Any]) -> pd.Series:
    start_depth = 0
    default_depth = 1000
    return_code = {'layers': 1, 'depths': [default_depth]}
    depth_builder = []
    depths_fillup = np.zeros(shape=(3,4), dtype=float)
    # col 0 layer is existing
    # col 1 depth countdown
    # col 2 depth collecting through analysis
    nlayers = 0
    SOL_ZMX = 0
    SOL_Z1 = 0
    SOL_Z2 = 0
    SOL_Z3 = 0
    SOL_Z4 = 0
    try:
        if isinstance(grammar_loimis, dict) and grammar_loimis.get('count', 0) <= 0:
            # logger.info('problem?')
            # return return_code
            pd.Series((nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4))
        # logger.info(result)
        # logger.info(type(result))
        soilparts_count = grammar_loimis['count']
        # col 1 depth countdown
        depths_fillup[1,0] = float(default_depth)
        
        
        for y in range(0, soilparts_count):
            constituents_obj = grammar_loimis['soilparts'][y]
            constituents_count = constituents_obj['count']
            constituents_arr = constituents_obj['constituents']
            # one soil part can have multiple constituents e.g. r and ls, but usually only one depth
            
            if len(constituents_arr) > 0 and not has_only_no_info_constituents(constituents_arr):
                # col 0 layer is existing
                depths_fillup[0,y] = 1

                depth_for_soilpart: Dict = {}
                depth_found_for_soilpart = False
                range_val = False
                from_val = 0
                to_val = 0
                depth_avg = 0
                for x in range(0, constituents_count):
                    constitu = constituents_arr[x]
                    if constitu['depth']:
                        depth_for_soilpart = constitu['depth']
                        depth_found_for_soilpart = True
                if depth_found_for_soilpart:
                    try:
                        range_val = depth_for_soilpart['range']
                        # if 'deeper_than' in depth_for_soilpart:
                        #     deeper_than_val = depth_for_soilpart['deeper_than']
                        if 'from' in depth_for_soilpart:
                            from_val = int(depth_for_soilpart['from'])
                        if 'to' in depth_for_soilpart:
                            to_val = int(depth_for_soilpart['to'])
                        if range_val:
                            depth_avg = int((to_val + from_val)/2)
                        else:
                            depth_avg = to_val

                        # col 2 depth collecting through analysis
                        depths_fillup[2,y] = float(depth_avg)

                    except KeyError as ex:
                        tb = traceback.format_exc()
                        logger.info("keyerror at {}".format(ex))
                        logger.info(tb)
                else:
                    down_to_def = default_depth - start_depth
                    depth_avg = down_to_def

                thickness = start_depth+depth_avg - start_depth

                if thickness <= 0:
                    # logger.info('layers thickness of layer/soilpart {} is zero( thickness{},depth_avg{},start_depth{} )'
                    #       .format(y+1, thickness,depth_avg,start_depth))

                    return_code.update({'layers': int(soilparts_count-1)})


                else:
                    if y == 0:
                        return_code.update({'depths': []})
                    depth_builder.append(int(thickness))
                    start_depth = start_depth+depth_avg
        
        # col 0 layer is existing
        set_layers = np.count_nonzero(depths_fillup[0,:])
        # col 2 depth collecting through analysis
        set_depths = np.count_nonzero(depths_fillup[2,:])
        sum_set_depths = depths_fillup[2,:].sum()
        
        if sum_set_depths < default_depth and set_depths < set_layers:
            spread_avg = float((default_depth - sum_set_depths) / (set_layers - set_depths))
            for i in range(0, set_layers):
                if depths_fillup[2,i] <= 0:
                    depths_fillup[2,i]= np.rint(spread_avg)
        elif sum_set_depths >= default_depth and set_depths < set_layers:
            num_diff_layers = set_layers - set_depths
            spread_avg = sum_set_depths / 10
            for i in range(0, set_layers):
                if depths_fillup[2,i] <= 0:
                    depths_fillup[2,i] = np.rint(spread_avg)
        
        return_code.update({'layers': int(set_layers)})
        i=0
        for x in range(0, set_layers):
            i=i+depths_fillup[2,x]
            return_code['depths'].append(i)
        # return_code.update({'depths': [i+d for d in depths_fillup[2,:]]})
        # result = return_code
        nlayers = int(return_code['layers'])
        depths = return_code['depths']
        SOL_ZMX = depths_fillup[2,:].sum()
        SOL_Z1 = depths[0]
        if len(depths) > 1 and nlayers > 1:
            SOL_Z2 = depths[1]
        if len(depths) > 2 and nlayers > 2:
            SOL_Z3 = depths[2]
        if len(depths) > 3 and nlayers > 3:
            SOL_Z4 = depths[3]
        return pd.Series((nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4))
    except KeyError as ex:
        tb = traceback.format_exc()
        logger.info("keyerror at {}".format(ex))
        logger.info(tb)
        return pd.Series((nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4))
    except Exception as ex:
        logger.info(ex)
        return pd.Series((nlayers, SOL_ZMX, SOL_Z1, SOL_Z2, SOL_Z3, SOL_Z4))


def try_texture_rules(idx1: str, idx2: str) -> Any:
    """tries to extract the correct numerical value from the lookup dict,
    'l2': {'sand': 90, 'silt': 3, 'clay': 7, 'lxtype': 'S'}

    Arguments:
        idx1 {str} -- main texture type index
        idx2 {str} -- sub-index for fraction type access

    Returns:
        Any -- mostly the integer value extracted for the type index
    """
    try:
        return texture_rules[idx1][idx2]
    except KeyError as ex:
        logger.warn('texture rule doesnt exist for {} // texture_rules[{}][{}]'.format(ex, idx1, idx2))
        return 0


def search_main_loimis_params(grammar_loimis: Dict[str, Any]) -> Dict[str, Any]:
    """extracts the texture fine earth and skeleton fractions from the data structure

    Arguments:
        grammar_loimis {grammar_loimis: Dict[str, Any]} -- the grammar data structure

    Returns:
        Dict[str, Any] -- a dict with fine earth and skeleton texture fields
    """

    peenes_builder: Dict = {
        'in_layers': []
    }

    kores_builder: Dict = {
        'in_layers': []
    }
        
    kores_types_builder: Dict = {
        'in_layers': []
    }

    return_code = {
        'layers': 0,
        'message': 'no_info'
    }

    try:
        if isinstance(grammar_loimis, dict) and grammar_loimis.get('count', 0) <= 0:
            # logger.info('problem?')
            return_code.update({'message': grammar_loimis.get('message', 'no_info')})
            grammar_loimis.get('count', 0)
            # peenes_builder.update({'1': 'no_peenes'})
            # peenes_builder.update({'in_layers': [1]})
            return_code.update({'layers': grammar_loimis.get('count', 0)})
            return_code.update({'kores': kores_builder})
            return_code.update({'kores_types': kores_builder})
            return_code.update({'peenes': peenes_builder})
            # logger.warn(return_code)
            return return_code

        result = grammar_loimis
        # logger.info(result)
        soilparts_count = grammar_loimis.get('count', 0)
        return_code.update({'layers': int(soilparts_count)})

        found_kores_overall = False
        found_main_peenes_overall = False

        for y in range(0, soilparts_count):

            found_kores = False
            found_kores_constitu_pos = -1
            found_main_peenes = False

            constituents_obj = result['soilparts'][y]
            constituents_count = constituents_obj['count']
            constituents_arr = constituents_obj['constituents']

            for x in range(0, constituents_count):
                constitu = constituents_arr[x]
                logger.debug(f"constitu: {constitu}")
                
                # find kores
                if constitu['type'] == 'kores' and found_kores == False:
                    
                    kores_type = constitu['code']
                    found_kores = True
                    found_kores_constitu_pos = x
                    found_kores_overall = True
                    
                    # logger.info(f"kores type: {kores_type}")
                    
                    kores_builder['in_layers'].append(y + 1)
                    kores_types_builder['in_layers'].append(y + 1)

                    if 'amp' in constitu and constitu['amp']:
                        kores_amp = constitu['amp']
                        
                        kores_builder.update({str(y + 1): kores_amp_strength[kores_amp]})
                        kores_types_builder.update({str(y + 1): kores_type})
                        # logger.info(f"kores_amp: {kores_amp}")

                    elif kores_type in ['pk', 'kr', 'p', 'd', 'lu', 'ck']:
                        # kr, p, d, lu, pk - astmeid ei kasutata
                        kores_builder.update({str(y + 1): int(kores_amp_strength[6])})
                        kores_types_builder.update({str(y + 1): kores_type})

                    elif kores_type in ['r⁰', 'r°', 'r' ,'v⁰', 'v°', 'v','kb','kb⁰','kb°','k⁰', 'k°', 'k']:
                        # mulla legend says if no amp, then highest kores
                        kores_builder.update({str(y + 1): int(kores_amp_strength[6])})
                        kores_types_builder.update({str(y + 1): kores_type})
                    else:
                        kores_builder.update({str(y + 1): 0})
                        kores_types_builder.update({str(y + 1): kores_type})

                # find peenes or turfs
                if (constitu['type'] == 'peenes' or constitu['type'] == 'turfs') and found_main_peenes == False:
                    peenes_type = constitu['code']
                    found_main_peenes = True
                    found_main_peenes_overall = True
                    
                    peenes_builder['in_layers'].append(y + 1)

                    if constitu['amp']:
                        main_peenes_amp = constitu['amp']
                        peenes_builder.update({str(y + 1): peenes_type + str(main_peenes_amp)})
                    else:
                        peenes_builder.update({str(y + 1): peenes_type})

            # new thing here to try capture no peenes but kores, means karst (or generally purely rocky)
            # this is still per layer
            if found_main_peenes == False and found_kores == True and found_kores_constitu_pos > -1:
                x = found_kores_constitu_pos
                constitu_x = constituents_arr[x]
                kores_type = constitu_x['code']

                peenes_builder['in_layers'].append(y + 1)

                if 'amp' in constitu_x and constitu_x['amp']:
                    kores_amp = constitu_x['amp']
                    # kores_builder.update({str(y + 1): kores_amp_strength[kores_amp]})
                    peenes_builder.update({str(y + 1): kores_type + str(kores_amp)})
                elif kores_type in ['no_info', 'pk', 'kr', 'p', 'd', 'lu', 'ck']:
                    # kr, p, d, lu, pk - astmeid ei kasutata
                    # kores_builder.update({str(y + 1): int(kores_amp_strength[3])})
                    peenes_builder.update({str(y + 1): kores_type})
                else:
                    # mulla legend says if no amp, then highest kores
                    # kores_builder.update({str(y + 1): int(kores_amp_strength[6])})
                    peenes_builder.update({str(y + 1): kores_type})

        message = 'kores_overall={}, peenes_overall={}'.format(
            found_kores_overall, found_main_peenes_overall)
        return_code.update({'message': message})
        return_code.update({'kores': kores_builder})
        return_code.update({'kores_types': kores_types_builder})
        return_code.update({'peenes': peenes_builder})
        logger.info(return_code)
        return return_code

    except KeyError as ex:
        tb = traceback.format_exc()
        logger.info("keyerror at {}".format(ex))
        logger.info(tb)
        return_code.update({'message': tb})
        peenes_builder.update({'1': 'no_info'})
        peenes_builder.update({'in_layers': [1]})
        return_code.update({'peenes': peenes_builder})
        return return_code
    except Exception as ex:
        logger.info('Exception: {}'.format(ex))
        return_code.update({'message': ex})
        peenes_builder.update({'1': 'no_info'})
        peenes_builder.update({'in_layers': [1]})
        return_code.update({'peenes': peenes_builder})
        return return_code


def set_texture_values_dask(grammar_loimis: Dict[str, Any]) -> pd.Series:
    """applies the extraction of the texture fine earth and skeleton fractions from the data structure to DataFrame row

    Arguments:
        grammar_loimis {grammar_loimis: Dict[str, Any]} -- the grammar data structure

    Returns:
        pd.Series -- a Pandas Series with fine earth and skeleton texture fields to be added to soil map DataFrame
    """

    result: Dict[str, Any] = search_main_loimis_params(grammar_loimis)
    num_layers = grammar_loimis.get('count', 0)
    
    
    
    actual_code = ''
    if isinstance(grammar_loimis, dict) and 'code' in grammar_loimis:
        actual_code = grammar_loimis.get('code','no_info')

    # (layer): Rock fragment content (% total weight) required
    sol_rock_container = [0, 0, 0, 0]
    rock_type_container = ['', '', '', '']

    if 'kores' in result and 'in_layers' in result['kores'] and len(result['kores']['in_layers']) > 0:
        # korrr = result['kores']['in_layers']
        # layer_num = korrr[0]
        korrr: List = result['kores'].get('in_layers')
        layer_num: str = str(korrr[0])
        
        sol_rock_container[int(layer_num) - 1] = result['kores'][layer_num]
        rock_type_container[int(layer_num) - 1] = result['kores_types'][layer_num]

        
        if len(korrr) > 1 and num_layers >= 1:
            layer_num: str = str(korrr[1])
            if int(layer_num) - 1 < 0 or int(layer_num) - 1 > 2:
                logger.error(
                    'kores layer_num out of range {}'.format(layer_num))
                logger.debug(
                    'result[kores]: {} - str(korrr[1]): {}'.format(result['kores'], korrr[1]))
                logger.debug(result)
            sol_rock_container[int(layer_num) - 1] = result['kores'][layer_num]
            rock_type_container[int(layer_num) - 1] = result['kores_types'][layer_num]
            
        if len(korrr) > 2 and num_layers >= 2:
            layer_num: str = str(korrr[2])
            sol_rock_container[int(layer_num) - 1] = result['kores'][layer_num]
            rock_type_container[int(layer_num) - 1] = result['kores_types'][layer_num]
            
        if len(korrr) > 3 and num_layers >= 3:
            layer_num: str = str(korrr[3])
            sol_rock_container[int(layer_num) - 1] = result['kores'][layer_num]
            rock_type_container[int(layer_num) - 1] = result['kores_types'][layer_num]
            

    SOL_ROCK1 = sol_rock_container[0]
    SOL_ROCK2 = sol_rock_container[1]
    SOL_ROCK3 = sol_rock_container[2]
    SOL_ROCK4 = sol_rock_container[3]
    
    EST_CRS1 = rock_type_container[0]
    EST_CRS2 = rock_type_container[1]
    EST_CRS3 = rock_type_container[2]
    EST_CRS4 = rock_type_container[3]
    
    SOL_CLAY1 = 0  # (layer ): Clay content (% soil weight) required
    SOL_CLAY2 = 0
    SOL_CLAY3 = 0
    SOL_CLAY4 = 0
    SOL_SILT1 = 0  # (layer ): Silt content (% soil weight) required
    SOL_SILT2 = 0
    SOL_SILT3 = 0
    SOL_SILT4 = 0
    SOL_SAND1 = 0  # (layer ): Sand content (%s soil weight) required
    SOL_SAND2 = 0
    SOL_SAND3 = 0
    SOL_SAND4 = 0
    LXTYPE1 = ''
    LXTYPE2 = ''
    LXTYPE3 = ''
    LXTYPE4 = ''

    EST_TXT1 = ''  # (layer ) only info to know which loimis per layer
    EST_TXT2 = ''
    EST_TXT3 = ''
    EST_TXT4 = ''
    
    if 'peenes' in result and 'in_layers' in result['peenes'] and len(result['peenes']['in_layers']) > 0:
        peen: List = result['peenes'].get('in_layers')
        layer_num: str = str(peen[0])
        idx = result['peenes'][layer_num]
        SOL_CLAY1 = try_texture_rules(idx, 'clay')
        SOL_SILT1 = try_texture_rules(idx, 'silt')
        SOL_SAND1 = try_texture_rules(idx, 'sand')
        LXTYPE1 = try_texture_rules(idx, 'lxtype')
        # maybe as additional info
        EST_TXT1 = idx
        if (SOL_CLAY1 + SOL_SILT1 + SOL_SAND1) != 100:
            logger.error('rules wrong for 1: {},{},{}'.format(actual_code, (SOL_CLAY1 + SOL_SILT1 + SOL_SAND1),
                                                             grammar_loimis['code']))

        if len(peen) > 1 and num_layers >= 1:
            idx = result['peenes'][str(peen[1])]
            SOL_CLAY2 = try_texture_rules(idx, 'clay')
            SOL_SILT2 = try_texture_rules(idx, 'silt')
            SOL_SAND2 = try_texture_rules(idx, 'sand')
            LXTYPE2 = try_texture_rules(idx, 'lxtype')
            EST_TXT2 = idx

            if (SOL_CLAY2 + SOL_SILT2 + SOL_SAND2) != 100:
                logger.error('rules wrong for 2: {},{},{}'.format(actual_code, (SOL_CLAY2 + SOL_SILT2 + SOL_SAND2),
                                                                 grammar_loimis['code']))
        if len(peen) > 2 and num_layers >= 2:
            idx = result['peenes'][str(peen[2])]
            SOL_CLAY3 = try_texture_rules(idx, 'clay')
            SOL_SILT3 = try_texture_rules(idx, 'silt')
            SOL_SAND3 = try_texture_rules(idx, 'sand')
            LXTYPE3 = try_texture_rules(idx, 'lxtype')
            EST_TXT3 = idx

            if (SOL_CLAY3 + SOL_SILT3 + SOL_SAND3) != 100:
                logger.info('rules  wrong for 3: {} | {} | {}'.format(actual_code, (SOL_CLAY3 + SOL_SILT3 + SOL_SAND3), grammar_loimis['code']))
        
        if len(peen) > 3 and num_layers >= 3:
            idx = result['peenes'][str(peen[3])]
            SOL_CLAY4 = try_texture_rules(idx, 'clay')
            SOL_SILT3 = try_texture_rules(idx, 'silt')
            SOL_SAND4 = try_texture_rules(idx, 'sand')
            LXTYPE4 = try_texture_rules(idx, 'lxtype')
            EST_TXT4 = idx

            if (SOL_CLAY4 + SOL_SILT4 + SOL_SAND4) != 100:
                logger.info('rules  wrong for 4: {} | {} | {}'.format(actual_code, (SOL_CLAY4 + SOL_SILT4 + SOL_SAND4), grammar_loimis['code']))

    # Texture field descript/code of soil layer optional (lets make siffer and main_EST_TXT?)
    # TEXTURE = full_texture
    # SNAM = row['simplified_siffer'] + '-' + EST_TXT1  # Soil Name
    # HYDGRP = set_hydrogrp(LXTYPE1)
    # we do that later in generllization,not needed yet
    # return TEXTURE, SNAM, HYDGRP

    return pd.Series((EST_TXT1, EST_CRS1, LXTYPE1, SOL_CLAY1, SOL_SILT1, SOL_SAND1, SOL_ROCK1,
                      EST_TXT2, EST_CRS2, LXTYPE2, SOL_CLAY2, SOL_SILT2, SOL_SAND2, SOL_ROCK2,
                      EST_TXT3, EST_CRS3, LXTYPE3, SOL_CLAY3, SOL_SILT3, SOL_SAND3, SOL_ROCK3,
                      EST_TXT4, EST_CRS4, LXTYPE4, SOL_CLAY4, SOL_SILT4, SOL_SAND4, SOL_ROCK4, result))


######################################################################
# HUMUS
###############################################3


class HumusVisitor(PTNodeVisitor):
    
    
    # def th_turfs_grammar(): return th_turfs, EOF
    # def t_turfs_grammar(): return t_turfs, EOF
    # def huumus_hor_grammar(): return huumus_hor, EOF
    # def kodu_depth_range_grammar(): return kodu_depth_range, EOF

    def visit_th_turfs_grammar(self, node, children):
        """
        th_turfs, EOF
        """
        if self.debug:
            logger.debug("th_turfs_grammar node {}".format(node.value))
            logger.debug("th_turfs_grammar {}".format(children))
        
        tmp = [p for p in children if isinstance(p, dict)]

        # we just reduce to one
        loimis_def = {"type": "th_turfs",
                      "count": len(tmp),
                      "soilparts": tmp}
        # return loimis_def
        return tmp[0]
    
    def visit_t_turfs_grammar(self, node, children):
        """
        t_turfs, EOF
        """
        if self.debug:
            logger.debug("visit_t_turfs_grammar node {}".format(node.value))
            logger.debug("visit_t_turfs_grammar {}".format(children))
            
        tmp = [p for p in children if isinstance(p, dict)]

        # we just reduce to one
        loimis_def = {"type": "t_turfs",
                      "count": len(tmp),
                      "soilparts": tmp}
        # return loimis_def
        return tmp[0]
    
    
    def visit_huumus_hor_grammar(self, node, children):
        """
        huumus_hor, EOF
        """
        if self.debug:
            logger.debug("huumus_hor_grammar node {}".format(node.value))
            logger.debug("huumus_hor_grammar {}".format(children))
            
        tmp = [p for p in children if isinstance(p, dict)]

        # we just reduce to one
        loimis_def = {"type": "huumus_hor",
                      "count": len(tmp),
                      "soilparts": tmp}
        # return loimis_def
        return tmp[0]
        
        
    def visit_kodu_depth_range_grammar(self, node, children):
        """
        kodu_depth_range, EOF
        """
        if self.debug:
            logger.debug("kodu_depth_range_grammar node {}".format(node.value))
            logger.debug("kodu_depth_range_grammar {}".format(children))
            
        tmp = [p for p in children if isinstance(p, dict)]

        # we just reduce to one
        loimis_def = {"type": "kodu_depth_range",
                      "count": len(tmp),
                      "soilparts": tmp}
        return loimis_def
            
            
    def visit_turfs(self, node, children):
        if self.debug:
            logger.debug("turfs node {}".format(node.value))
            logger.debug("turfs {}".format(children))

        # contains turfs variants incl amps, depth
        dat = {"type": "turfs"}
        t_child = children
        
        # dat.update({"code": node.value})

        # the compulsory peenes/turvas code
        dat.update({"type": children[0]})
        t_child.pop(0)

        # an optional amp for the peenes
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'amp' in t_child[0]:
            dat.update(t_child[0])
            t_child.pop(0)
        else:
            dat.update({"amp": False})

        # an optional depth related
        if len(t_child) > 0 and isinstance(t_child[0], dict) == True and 'depth' in t_child[0]:
            dat.update(t_child[0])
        else:
            dat.update({"depth": False})

        return dat

    
    def visit_th_turfs(self, node, children):
        """
        th, depth_range
        """
        return self.visit_turfs(node, children)
    
    def visit_t_turfs(self, node, children):
        """
        t, Optional(amplifiers), depth_range
        """
        return self.visit_turfs(node, children)
    

    
    def visit_kodu_depth_range(self, node, children):
        """
        kodu_depth_number, ZeroOrMore('-', kodu_depth_number)
        """
        if self.debug:
            logger.debug("visit_kodu_depth_range node {}".format(node.value))
            logger.debug("visit_kodu_depth_range {}".format(children))
        
        dat = {"type": "kodu",
              "depth": {}
              }
        
        t_child = children
        
        if len(t_child) == 2:
            logger.info("depth is from to {}".format(children))
            dat['depth'].update(
                {"range": True,
                 "from": t_child[0],
                 "to": t_child[1]
                 })

        if len(t_child) == 1:
            dat['depth'].update(
                {"range": False,
                 "to": t_child[0]
                 })
        # compound operator for depth elements, either has one number, or two numbers separated with '-'
        # should return the calculated depth horizon?
        return dat
        

    def visit_huumus_hor(self, node, children):
        """
        depth_range
        """
        if self.debug:
            logger.debug("visit_huumus_hor node {}".format(node.value))
            logger.debug("visit_huumus_hor {}".format(children))
                # contains turfs variants incl amps, depth
                
        dat = {"type": "plain_huumus",
              'depth': {}}
        
        t_child = children
        
        if isinstance(t_child, list):
            if len(t_child) == 2:
                dat['depth'].update(
                    {"range": True,
                     "from": t_child[0],
                     "to": t_child[1]
                     })

            if len(t_child) == 1:
                dat['depth'].update(
                    {"range": False,
                     "to": t_child[0]
                     })
        elif isinstance(t_child, dict):
            dat.update(t_child)
        else:
            dat.update({"depth": False})

        return dat
    
        
    def visit_t(self, node, children):
        return node.value
    
    
    def visit_th(self, node, children):
        return node.value
    
    
    def visit_amplifiers(self, node, children):
        """
        def amplifiers(): return [amp1, amp2, amp3, amp4, amp5]
        """
        if self.debug:
            logger.debug("amplifiers node {}".format(node.value))
            logger.debug("amplifiers {}".format(children))

        # compound list of amps, should return only one number 1-5, might not be needed for texture or similar?
        if node.value == '₁':
            return {"amp": 1}
        if node.value == '₂':
            return {"amp": 2}
        if node.value == '₃':
            return {"amp": 3}
        if node.value == '₄':
            return {"amp": 4}
        if node.value == '₅':
            return {"amp": 5}
        return {"amp": int(node.value)}

        
    def visit_kodu_depth_number(self, node, children):
        """
        depth_number, Optional(amplifiers)
        """
        if self.debug:
            logger.debug("kodu_depth_number {}".format(children))
            logger.debug("kodu_depth_number node {}".format(node.value))

        dat = {"depth": {}}
        t_child = children

        if len(t_child) == 2:
            # logger.info("depth has amp {}".format(children))
            dat['depth'].update(
                {"to": t_child[0],
                 "amp": t_child[1]
                 })

        if len(t_child) == 1:
            dat['depth'].update(
                {"to": t_child[0]
                 })
            
        return dat 
    
    
    def visit_depth_range(self, node, children):
        """
        depth_number, ZeroOrMore('-', depth_number)
        """
        if self.debug:
            logger.debug("depth_range {}".format(children))
            logger.debug("depth_range node {}".format(node.value))

        dat = {"depth": {}}
        t_child = children

        if len(t_child) == 2:
            # logger.info("depth is from to {}".format(children))
            dat['depth'].update(
                {"range": True,
                 "from": t_child[0],
                 "to": t_child[1]
                 })

        if len(t_child) == 1:
            dat['depth'].update(
                {"range": False,
                 "to": t_child[0]
                 })
        # compound operator for depth elements, either has one number, or two numbers separated with '-'
        # should return the calculated depth horizon?
        return dat


    def visit_depth_number(self, node, children):
        """
        def depth_number(): return RegExMatch(r'\\d\\d+') # at least two digits for depth numbers
        """
        if self.debug:
            logger.debug("depth_number {}".format(children))
        # itself might not be needed if handled already by range?
        # here we do the SWAT cm to mm conversion
        return float(node.value)*10


def get_humus_horizon(parse_tree):
    result_parsed = visit_parse_tree(parse_tree, HumusVisitor(debug=False))
    # pprint.pprint(result_parsed)
    depth_value_avg = 0.0
    # pprint.pprint(result_parsed.get('depth',{}))
    depths = result_parsed.get('depth',{})
    if depths.get('range', False) is False:
        depth_value_avg = depths.get('to', 0.0)
    else:
        from_v = depths.get('from', 0.0)
        to_v = depths.get('to', 0.0)
        if to_v > 0 and from_v > 0 and to_v > from_v:
            depth_value_avg = np.rint( (to_v + from_v) /  2 )
        else:
            logger.debug(result_parsed)
            # could be just 0
            
    horizon_type = result_parsed.get('type','no_info')
    return (depth_value_avg, horizon_type, result_parsed)


def analyse_huumus(x, parsers_d):
    if len(x.split("/")) > 1:
        lower_horizon = x.split("/")[1]
        p = can_parse_multiple_get_parser(lower_horizon, parsers_d)
        if p is not None:
            t = p.parse(lower_horizon)
            # logger.info("has two layers: {} >> parse_tree: {}".format(lower_horizon, t))
            gr = get_humus_horizon(t)
            return pd.Series((gr[0], gr[1], 'lower_h', gr[2]))
        else:
            # depth_value_avg, horizon_type, 'single_h', result_parsed
            return pd.Series((0, 'no_info', 'single_h', {'error': x}))
    else:
        sinlge_horizon = x
        p = can_parse_multiple_get_parser(sinlge_horizon, parsers_d)
        if p is not None:
            t = p.parse(sinlge_horizon)
            # logger.info("single_horizon: {} >> parse_tree: {}".format(x, t))
            gr = get_humus_horizon(t)
            return pd.Series((gr[0], gr[1], 'single_h', gr[2]))
        else:
            # depth_value_avg, horizon_type, 'single_h', result_parsed
            return pd.Series((0, 'no_info', 'single_h', {'error': x}))


def sort_ao_horizon_type(row, fieldname_select, fieldname_return, selecttype):
    retval = 0.0
    if row[fieldname_select] == selecttype:
        return row[fieldname_return]
    else:
        return retval



def demo_soc_stock_i_MgC_ha_for_layer(row, layernum):
    if int(row['nlayers']) <  layernum:
        return pd.Series((0, 0, 0))
    depth_i = row[f'SOL_Z{layernum}']
    # m2_area = row['area']
    # ha_area=m2_area/10000
    # cm2_area = m2_area*10000
    if layernum > 1:
        for t_num in range(1, layernum):
            depth_i = depth_i - row[f'SOL_Z{t_num}']
    
    # depth is in mm
    # volume here is in cm
    cm_depth = depth_i/10
    
    BDfinei = row[f'sol_bd{layernum}']
    
    BD_manual = ((row[f'SOL_CLAY{layernum}']/100*clay_dg+row[f'SOL_SILT{layernum}']/100*loam_dg+row[f'SOL_SAND{layernum}']/100*sand_dg) - row[f'SOL_ROCK{layernum}']/100*loam_dg) / (1 - row[f'SOL_ROCK{layernum}']/100)
    
    soc_stock_i_MgC_ha = row[f'sol_soc{layernum}'] * BDfinei * (1 - row[f'SOL_ROCK{layernum}']/100 ) * cm_depth * 0.1 * ha_area
    soc_stock_i_MgC_ha_manual = row[f'sol_soc{layernum}'] * BD_manual * (1 - row[f'SOL_ROCK{layernum}']/100 ) * cm_depth * 0.1 * ha_area
    
    return pd.Series((BD_manual, soc_stock_i_MgC_ha, soc_stock_i_MgC_ha_manual))


def soc_stock_i_MgC_ha_for_layer(row, layernum):
    if int(row['nlayers']) <  layernum:
        return 0.0
    depth_i = row[f'SOL_Z{layernum}']
    
    if layernum > 1:
        depth_i = depth_i - row[f'SOL_Z{layernum-1}']
    
    # depth is in mm
    # volume here is in cm
    cm_depth = depth_i/10
    
    BDfinei = row[f'sol_bd{layernum}']
    soc_stock_i_MgC_ha = float(row[f'sol_soc{layernum}'] * BDfinei * (1 - row[f'SOL_ROCK{layernum}']/100 ) * cm_depth * 0.1)
    return soc_stock_i_MgC_ha


def avg_exisiting_layers(row):
    l = int(row['nlayers'])
    if l == 1:
        return row['soc_stock_MgC_ha_1']
    elif l == 2:
        return (row['soc_stock_MgC_ha_1'] + row['soc_stock_MgC_ha_2'] ) / 2
    elif l == 3:
        return (row['soc_stock_MgC_ha_1'] + row['soc_stock_MgC_ha_2'] + row['soc_stock_MgC_ha_3']) / 3
    elif l == 4:
        return (row['soc_stock_MgC_ha_1'] + row['soc_stock_MgC_ha_2'] + row['soc_stock_MgC_ha_3'] + row['soc_stock_MgC_ha_4']) / 4
    else:
        return 0.0


