"""psOperatorSpecializer.py -- tool to specialize genearal T2 operators"""

import collections
import os
import sys

from fontTools.ttLib import TTFont

def specializeOperators(font):
    cff = font['CFF '].cff

    for td in cff.topDictIndex:
        for cs in td.CharStrings.values():
            cs.decompile()
            cs.program = simplifyCharstring(cs.program)
        for cs in td.GlobalSubrs:
            if cs.program:
                cs.program = simplifyCharstring(cs.program)
        if hasattr(td, 'FDArray'):
            for fd in td.FDArray:
                if hasattr(fd, 'Subrs'):
                    for cs in fd.Subrs:
                        if cs.program:
                            cs.program = simplifyCharstring(cs.program)

def group(el, n):
    """Group the list el into groups of size n"""

    for i in xrange(0, len(el), n):
        yield el[i:i+n]

def formatter_hhcurveto(dxa, dya, dxb, dyb, dxc, dyc):
    if dya != 0:
        return [dya, dxa, dxb, dyb, dxc]
    else:
        return [dxa, dxb, dyb, dxc]

def formatter_hvcurveto(dxa, dya, dxb, dyb, dxc, dyc):
    if dya == 0:
        # second "should" be vertical
        if dxc != 0:
            return [dxa, dxb, dyb, dyc, dxc]
        else:
            return [dxa, dxb, dyb, dyc]
    else:
        # second "should" be horizontal
        if dyc != 0:
            return [dya, dxb, dyb, dxc, dyc]
        else:
            return [dya, dxb, dyb, dxc] 

formatter_vhcurveto = formatter_hvcurveto

def formatter_vvcurveto(dxa, dya, dxb, dyb, dxc, dyc):
    if dxa != 0:
        return [dxa, dya, dxb, dyb, dyc]
    else:
        return [dya, dxb, dyb, dyc]

def simplifyCharstring(old_program):
    new_program = []
    stack = []

    for tok in old_program:
        if not isinstance(tok, basestring):
            stack.append(tok)
        elif tok == 'rmoveto':
            if len(stack) == 3 and not new_program:
                # this is a width arg
                new_program.append(stack[0])
                del stack[0]
            else:
                assert len(stack) == 2
            if stack[0] == 0 and stack[1] == 0:
                # it's a no-op
                pass
            elif stack[0] == 0 and stack[1] != 0:
                new_program.extend([stack[1], 'vmoveto'])
            elif stack[1] == 0:
                new_program.extend([stack[0], 'hmoveto'])
            else:
                new_program.extend(stack + [tok])
            stack = []
        elif tok == 'rlineto':
            assert len(stack) % 2 == 0 and len(stack) >= 2

            if len(stack) == 2 and new_program[-1] == 'rrcurveto':
                new_program[-1:] = stack + ['rcurveline']
                stack = []
                continue

            # get (x, y) pairs
            groupgen = group(stack, 2)
            toksets = [(x, y) for x, y in groupgen if x != 0 or y != 0]

            # join back-to-back horizontal or vertical lines
            new_toksets = toksets[:1]
            for cur in toksets[1:]:
                for i, other in zip(range(2), reversed(range(2))):
                    if new_toksets[-1][i] == 0 and cur[i] == 0:
                        new_toksets[-1][other] += cur[other]
                    else:
                        new_toksets.append(cur)

            # split into ranges
            ranges = []
            start = 0
            is_hv = toksets[0][0] == 0 or toksets[0][1] == 0
            tsiter = iter(enumerate(toksets))
            next(tsiter)
            for idx, (x, y) in tsiter:
                if not is_hv and (x == 0 or y == 0):
                    ranges.append((start, idx, False))
                    start = idx
                    is_hv = True
                elif is_hv and (x != 0 and y != 0):
                    ranges.append((start, idx, True))
                    start = idx
                    is_hv = False
            ranges.append((start, len(toksets), is_hv))
            if ranges[0][1] - ranges[0][0] == 0:
                del ranges[0]

            # process into output
            needs_close = False
            for start, stop, is_hv in ranges:
                if is_hv and ((not needs_close and stop - start > 1) or (needs_close and stop - start > 2)):
                    if needs_close:
                        new_program.append('rlineto')
                    on_x = toksets[start][0] != 0
                    for x, y in toksets[start:stop]:
                        if on_x:
                            new_program.append(x)
                        else:
                            new_program.append(y)
                        on_x = not on_x
                    if toksets[start][0] == 0:
                        new_program.append('vlineto')
                    else:
                        new_program.append('hlineto')
                else:
                    new_program.extend(stack[start * 2:stop * 2])
                    needs_close = True
            if not isinstance(new_program[-1], basestring):
                new_program.append('rlineto')

            stack = []
        # TODO(sfishman): hlineto and vlineto
        elif tok == 'rrcurveto':
            assert len(stack) % 6 == 0 and len(stack) >= 6

            # group into curves
            groupgen = group(stack, 6)
            toksets = list(groupgen)

            # dp to figure out best encoding
            results = [0 for _ in xrange(len(toksets))]
            split = [0 for _ in xrange(len(toksets))]
            op = [0 for _ in xrange(len(toksets))]
            for i in xrange(len(toksets)):
                # encode last item seperately from the rest
                min_option = 7 + results[i - 1]
                min_k = i
                min_oper = 'rrcurveto'
                for k in range(i + 1):
                    # consider grouping k..i (inclusive)
                    option = (i - k + 1) * 6 + 1 # rrcurveto
                    poss = {'hhcurveto': True, 'hvcurveto': True,
                            'vhcurveto': True, 'vvcurveto': True}
                    if toksets[k][5] != 0:
                        poss['hhcurveto'] = False
                    if toksets[k][1] != 0:
                        poss['hvcurveto'] = False
                    if toksets[k][0] != 0:
                        poss['vhcurveto'] = False
                    if toksets[k][4] != 0:
                        poss['vvcurveto'] = False
                    alt = True
                    for j in xrange(k + 1, i):
                        if toksets[j][1] != 0:
                            poss['hhcurveto'] = False
                            if not alt:
                                poss['hvcurveto'] = False
                            else:
                                poss['vhcurveto'] = False
                        if toksets[j][0] != 0:
                            poss['vvcurveto'] = False
                            if not alt:
                                poss['vhcurveto'] = False
                            else:
                                poss['hvcurveto'] = False
                        if toksets[j][4] != 0:
                            poss['vvcurveto'] = False
                            if not alt:
                                poss['hvcurveto'] = False
                            else:
                                poss['vhcurveto'] = False
                        if toksets[j][5] != 0:
                            poss['hhcurveto'] = False
                            if not alt:
                                poss['vhcurveto'] = False
                            else:
                                poss['hvcurveto'] = False
                        if not all(poss):
                            break
                    if i != k:
                        if toksets[i][5] != 0:
                            poss['hhcurveto'] = False
                        if toksets[i][1] != 0:
                            if not alt:
                                poss['hvcurveto'] = False
                            else:
                                poss['vhcurveto'] = False
                            poss['hhcurveto'] = False
                        if toksets[i][0] != 0:
                            if not alt:
                                poss['vhcurveto'] = False
                            else:
                                poss['hvcurveto'] = False
                            poss['vvcurveto'] = False
                        if toksets[i][4] != 0:
                            poss['vvcurveto'] = False

                    oper = 'rrcurveto'
                    if poss['hhcurveto']:
                        val = (i - k + 1) * 4 + 1
                        if toksets[k][1] != 0:
                            val += 1
                        if val < option:
                            oper = 'hhcurveto'
                            option = val
                    if poss['hvcurveto']:
                        val = (i - k + 1) * 4 + 1
                        if (i - k + 1) % 2 == 0:
                            # should end horizontal
                            if toksets[i][5] != 0:
                                val += 1
                        else:
                            # should end vertical
                            if toksets[i][4] != 0:
                                val += 1
                        if val < option:
                            oper = 'hvcurveto'
                            option = val
                    if poss['vhcurveto']:
                        val = (i - k + 1) * 4 + 1
                        if (i - k + 1) % 2 == 0:
                            # should end vertical
                            if toksets[i][4] != 0:
                                val += 1
                        else:
                            # should end horizontal
                            if toksets[i][5] != 0:
                                val += 1
                        if val < option:
                            oper = 'vhcurveto'
                            option = val
                    if poss['vvcurveto']:
                        val = (i - k + 1) * 4 + 1
                        if toksets[k][0] != 0:
                            val += 1
                        if val < option:
                            oper = 'vvcurveto'
                            option = val

                    option += results[k - 1]

                    if option < min_option:
                        min_option = option
                        min_k = k
                        min_oper = oper

                results[i] = min_option
                split[i] = min_k
                op[i] = min_oper

            # apply the results
            reform = collections.deque()
            cur_end = len(toksets) - 1
            while cur_end >= 0:
                cur_start = split[cur_end]
                cur_op = op[cur_end]

                segment = collections.deque()
                if cur_op == 'rrcurveto':
                    formatter = lambda *x: x
                elif cur_op == 'hhcurveto':
                    formatter = formatter_hhcurveto
                elif cur_op == 'hvcurveto':
                    formatter = formatter_hvcurveto
                elif cur_op == 'vhcurveto':
                    formatter = formatter_vhcurveto
                elif cur_op == 'vvcurveto':
                    formatter = formatter_vvcurveto
                for i in xrange(cur_start, cur_end + 1):
                    it = formatter(*toksets[i])
                    segment.extendleft(it)
                segment.appendleft(cur_op)
                reform.extendleft(segment)

                cur_end = cur_start - 1

            new_program.extend(reform)
            stack = []
        else:
            new_program.extend(stack)
            new_program.append(tok)
            stack = []

    return new_program
