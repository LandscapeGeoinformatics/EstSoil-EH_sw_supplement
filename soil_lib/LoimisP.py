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

from arpeggio import (EOF, And, Combine, NoMatch, Not, OneOrMore, Optional,  # type: ignore
                      PTNodeVisitor, RegExMatch, UnorderedGroup, ZeroOrMore)  # type: ignore

import logging
from operator import itemgetter

# create logger
logger = logging.getLogger(__name__)


class LoimispVisitor(PTNodeVisitor):
    """The main Arpeggio Visitor implementaion to evaluate the found expressions and build a nice data structure
    
    Arguments:
        PTNodeVisitor -- the class to be inherited from Arpeggio

    Returns:
        LoimispVisitor
    """


    def cut_list_after(self, the_list, separator):
        """

        :param the_list:
        :param separator:
        :return:
        """
        
        tmp = the_list
        if len(the_list) > 1:
            try:
                idx = the_list.index(separator)
                tmp = the_list[:idx]
            except ValueError as ex:
                logger.debug(ex)
                return tmp
        return tmp


    def visit_loimisp(self, node, children):
        """
        main grammar, visiting ...

        def loimisp(): return OneOrMore(soilParts, sep=horiSep), EOF
        """
        if self.debug:
            logger.debug("Loimis node {}".format(node.value))
            logger.debug("Loimis {}".format(children))

        # attention, SWAT needs depths in mm, but mullakaard are reported in cm
        # default_max_depth = 1000
        # default_stoniness = 2
        # layer stuff
        # clay, silt, sand percentages
        # rusle
        # soil awailable water content
        # organic carbon?!
        # drainage class A-D
        # has soilparts, but only take first and drop after first horiSep
        # here split and remove for horiSep
        tmp = self.cut_list_after(children, ';')
        if len(tmp) > 1:
            logger.info("loimisp untypically long after horisep cut({}) {}".format(
                len(children), children))

        # we just reduce to one
        loimis_def = {"type": "loimis",
                      "count": len(tmp[0]),
                      "soilparts": tmp[0]}

        return loimis_def

    def visit_soilParts(self, node, children):
        """
        def soilParts(): return paraSeq, ZeroOrMore(vertiSep, paraSeq)
        """
        if self.debug:
            logger.debug("soilParts node {}".format(node.value))
            logger.debug("soilParts({}) {}".format(len(children), children))
        # soil parts is the sum list of the sequence of soil parts (paraseq),
        # separated by a vertiSep (at different depths)
        # here split if has vertiSep
        tmp = [p for p in children if isinstance(p, dict)]
        return tmp

    def visit_paraSeq(self, node, children):
        """
        def paraSeq(): return OneOrMore(constituent)
        """
        if self.debug:
            logger.debug("paraSeq node {}".format(node.value))
            logger.debug("paraSeq({}) {}".format(len(children), children))

        # paraSeq is a sequence of constituents
        # if kores and peenes/turfs are in wrong order, multiple paraseqs can occur, maybe reorder?
        tmp = children
        if len(children) > 1:
            logger.info("paraSeq untypically long({}) {}".format(
                len(children), children))
            # we could reduce to just the first element? -> if first paraseq is complete (kores and peenes/turfs)
            # Possibly check if that's just wrong order (first peenes and second kores?) and make it one paraseq
            paraseq_final = []
            order_final = []
            used_elems = []
            paraseq_idx = 0
            has_kores = False
            has_peenes = False

            for para in children:
                con_elems = para['constituents']

                count_kores = len(
                    [k for k in con_elems if k['type'] == 'kores'])
                count_peenes = len(
                    [k for k in con_elems if k['type'] == 'peenes'])
                count_turfs = len(
                    [k for k in con_elems if k['type'] == 'turfs'])
                order = []

                for x in para['constituents']:
                    if isinstance(x, dict) and x['type'] == 'kores' and has_kores == False:
                        has_kores = True
                        order.append('K' + str(paraseq_idx))
                        used_elems.append('K' + str(paraseq_idx))
                        paraseq_final.append(x)
                    elif isinstance(x, dict) and x['type'] == 'kores' and has_kores == True:
                        order.append('K' + str(paraseq_idx))

                    if isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == False:
                        has_peenes = True
                        order.append('P' + str(paraseq_idx))
                        used_elems.append('P' + str(paraseq_idx))
                        paraseq_final.append(x)
                    elif isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == True:
                        order.append('P' + str(paraseq_idx))

                logger.info("paraseq_idx({}) order: _{}_ count_kores({}) count_peenes({}) count_turfs({})".format(
                    paraseq_idx, '_'.join(order), count_kores, count_peenes, count_turfs))
                paraseq_idx = paraseq_idx + 1

            logger.info("reused elements: _{}_ ".format('_'.join(used_elems)))
            sorted_constituents = sorted(paraseq_final, key=itemgetter('type'))
            len_constituents = len(sorted_constituents)
            paraseq_upd = {"count": len_constituents,
                           "constituents": sorted_constituents}
            tmp = [paraseq_upd]

            # debug final
            order_final = []
            for x in tmp[0]["constituents"]:
                if isinstance(x, dict) and x['type'] == 'kores':
                    order_final.append('K')

                if isinstance(x, dict) and x['type'] == 'peenes':
                    order_final.append('P')

            logger.info("final order: _{}_ ".format('_'.join(order_final)))

        return {"count": len(tmp),
                "paraseq": tmp}

    def visit_constituent(self, node, children):
        """
        def constituent(): return ZeroOrMore(kores), ZeroOrMore([peenes, turfs, Optional(alternateComma)])
        """
        if self.debug:
            logger.debug("constituent node {}".format(node.value))
            logger.debug("constituent({}) {}".format(len(children), children))

            if len(children) > 2:
                logger.info("constituent untypically long({}) {}".format(
                    len(children), children))

        # one compound element i.e. one type combo of kores, turf,soil with amps or karbonated and depth,
        # can have alternate comma sep, should stuff drop behind comma
        tmp_comma = self.cut_list_after(children, ',')
        tmp = self.cut_list_after(tmp_comma, '-')

        if len(tmp) > 2:
            logger.info("constituent untypically long after comma cut({}) {}".format(
                len(children), children))
            # we could reduce to one skeleton (koores) and one peat (turvas)/fine earth (peenes)?
            # typically kores-kores-peenes, kores-peenes-peenes
            # basically check for repetations (+1 same?) and pop the following, keep first

            count_kores = len([k for k in tmp if k['type'] == 'kores'])
            count_peenes = len([k for k in tmp if k['type'] == 'peenes'])
            count_turfs = len([k for k in tmp if k['type'] == 'turfs'])

            order = []
            has_kores = False
            has_peenes = False
            constituents_final = []
            for x in tmp:
                if isinstance(x, dict) and x['type'] == 'kores' and has_kores == False:
                    has_kores = True
                    order.append('K')
                    constituents_final.append(x)
                elif isinstance(x, dict) and x['type'] == 'kores' and has_kores == True:
                    order.append('K')

                if isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == False:
                    has_peenes = True
                    order.append('P')
                    constituents_final.append(x)
                elif isinstance(x, dict) and x['type'] == 'peenes' and has_peenes == True:
                    order.append('P')

            logger.debug("order: _{}_ count_kores({}) count_peenes({}) count_turfs({})"
                         .format('_'.join(order), count_kores, count_peenes, count_turfs))

            tmp = constituents_final

            # debug final
            order_final = []
            for x in constituents_final:
                if isinstance(x, dict) and x['type'] == 'kores':
                    order_final.append('K')

                if isinstance(x, dict) and x['type'] == 'peenes':
                    order_final.append('P')

            logger.info("final order: _{}_ ".format('_'.join(order_final)))

        return {"count": len(tmp),
                "constituents": tmp}

    def visit_kores(self, node, children):
        """
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

        # the compulsory peenes code
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

    def visit_kores_list(self, node, children):
        """
        def kores_list(): return [kr, kb, pk, ck, lu, v_0, k_0, r, v, k, p, d]
        """
        if self.debug:
            logger.debug("kores_list node {}".format(node.value))
            logger.debug("kores_list {}".format(children))
        # contains single kores variants from list
        if len(children) == 1:
            return children[0]
        else:
            logger.info("warning too many kores_list here {}".format(children))
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
        def depth_range(): return depth_number, ZeroOrMore('-', depth_number)
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

    def visit_horiSep(self, node, children):
        if self.debug:
            logger.debug("horiSep {}".format(children))
        # the horisep semicolon indicates another soil osa column within same polygon, skipped
        # indicator should be used at appropriate higher place to skip rest of sequence
        return node.value

    def visit_alternateComma(self, node, children):
        if self.debug:
            logger.debug("alternateComma {}".format(children))
        # alternateComma, indicator should be used at appropriate higher place to skip rest of sequence
        return node.value


def loimis_main_grammar():
    """defines the main Arpeggio grammar
    
    Returns:
        Arpeggio Grammar (Repetion?) -- needs an Arpeggio ParserPython object to initiate
    """


    def l(): return 'l'  # liiv, en: sand

    def pl(): return 'pl'  # pl - peenliiv, en: fine sand (täiendina peenliivakas)

    def plsl(): return 'plsl'  # plsl - peenliivakas saviliiv, en: fine clayey sand

    def sl(): return 'sl'  # sl - saviliiv, en: clayey sand

    def tsl(): return 'tsl'  # tsl - tolmjas saviliiv, en: dusty clayey sand

    def dk(): return 'dk'  # dk - liivakivirähk, en: sandstone grit

    def ls(): return 'ls'

    # def l(): return '' # ls₁ - kerge liivsavi, en: light sandy clay
    # def l(): return '' # ls₂ - keskmine liivsavi, en: medium sandy clay
    # def l(): return '' # ls₃ - raske liivsavi, en: heavy sandy clay

    def tls(): return 'tls'  # tls - tolmjas liivsavi, en: dusty sandy clay

    def s(): return 's'  # s - savi, en: clay

    def peenes_list(): return [plsl, pl, tsl, tls, dk, sl, ls, s, l]

    def peenes(): return Optional(kPlus), peenes_list, Optional(
        amplifiers), Optional(depth_range)

    # th 15 või th 15-20 toorhuumusliku horisondi tüsedus, en: raw humus thickness
    def th(): return 'th'

    def t(): return 't'  # t - turvas, en: peat

    def turfs(): return Optional(kPlus), [th, t], Optional(
        amplifiers), Optional(depth_range)

    def kr(): return 'kr'  # kr - kruus, en: gravel

    def r_0(): return ['r⁰', 'r°']  # r - rähk, en: grit, rubble

    def r(): return 'r'  # r - rähk, en: grit, rubble r⁰

    def v_0(): return ['v⁰', 'v°']  # v° - raudkiviveeris, en: granite rock

    def v(): return 'v'  # v - paeveeris, en: limestone rock

    def kb(): return 'kb'  # kb - klibu, en: debris, rubble

    def ck(): return 'ck'  # ck - kiltkivirähk, en: debris, rubble of slate/schist

    def k_0(): return ['k⁰', 'k°']  # k⁰ - raudkivid, en: granite rock

    def k(): return 'k'  # k - paekivid, en: limestone rock

    def pk(): return 'pk'  # pk - paeplaadid, en: limestone rock plates?

    def p(): return 'p'  # p - paas, en: limestone rock

    def d(): return 'd'  # d - liivakivi, en: sandstone rock?

    # lu - lubisetted, calcic deposits? kores peenes turfs alternateComma amplifiers amp1 depth_range
    def lu(): return 'lu'

    def kores_list(): return [kr, kb, pk, ck, lu, v_0, k_0, r_0, r, v, k, p, d]

    def kores(): return Optional(kPlus), kores_list, Optional(
        amplifiers), Optional(depth_range)

    # "+" r, v, kb,k, pk, p, lu - karbonaatne(kihisev materjal), but seems also to be used in more contexts
    def kPlus(): return '+'

    # can only be one digit (should be subscript, otherwise error), diff to depth numbers
    def amp1(): return Combine(['₁', (Not(depth_range), '1')])

    def amp2(): return Combine(['₂', (Not(depth_range), '2')])

    def amp3(): return Combine(['₃', (Not(depth_range), '3')])

    def amp4(): return '₄'

    def amp5(): return '₅'

    def amplifiers(): return [amp1, amp2, amp3, amp4, amp5]

    # at least two digits for depth numbers (well might as well be one)
    def depth_number(): return RegExMatch(r'\d+')

    def depth_range(): return Optional(
        kPlus), depth_number, ZeroOrMore('-', depth_number)

    def vertiSep(): return '/'

    def horiSep(): return ';'

    def alternateComma(): return [',', '-']

    def constituent(): return ZeroOrMore(kores), ZeroOrMore(
        [peenes, turfs, Optional(alternateComma)])

    def paraSeq(): return OneOrMore(constituent)

    def soilParts(): return paraSeq, ZeroOrMore(vertiSep, paraSeq)

    def loimisp(): return OneOrMore(soilParts, sep=horiSep), EOF

    return loimisp
