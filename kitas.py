from copy import deepcopy
from math import fabs, log, inf, ceil
from itertools import product

class KITAS:
    def __init__(self, eos, modulo_range, dimension=None):
        self.eos = eos
        self.modulo_range = modulo_range
        self.dimension = dimension

    def is_goal(self, act, epsilon, min_zeros, elem):
        if self.eos is not False and act is not None and \
                (act[:3] == 'ANA' or 'SYM' in act or act[:4] == 'RECK'):
            return True
        try:
            if min_zeros == 1:
                return False if None in self.eos[-1] else \
                    fabs(self.eos[-1][-1]-elem) < epsilon
            for i in self.eos[-min_zeros:]:
                for j in i:
                    if j is None or fabs(j-elem) >= epsilon:
                        return False
            return True
        except:
            return False

    def check_number(self, elem):
        if self.modulo_range == 10:
            return elem
        dec = elem-int(elem)
        elem = int(elem) if dec == 0.0 else elem
        if type(elem) == float:
            return False
        while elem >= self.modulo_range:
            elem -= self.modulo_range
        while elem < 0:
            elem += self.modulo_range
        return elem

    def do_op(self, op, m, n, rev):
        if op == '+' or (op == '-' and rev):
            return self.check_number(m + n)
        elif op == '-' or (op == '+' and rev):
            return self.check_number(m - n)
        elif op == '*' or (op == '/' and rev):
            return self.check_number(m * n)
        else:
            return self.check_number(m / n) if n != 0 else False

    def basic(self):
        if len(self.eos) < 2:
            return False
        t = [self.eos]
        for y in range(len(self.eos)-1):
            seq = []
            for x in range(len(t[y])-1):
                e = self.check_number(t[y][x+1]-t[y][x])
                if e is False:
                    return False
                seq.append(e)
            t.append(seq)
        return KITAS(t, self.modulo_range, self.dimension)

    def divisions(self):
        if len(self.eos[0]) < 2 or 0 in self.eos[0]:
            return False
        t = [self.eos[0]]
        for y in range(len(self.eos)-1):
            seq = []
            for x in range(len(t[y])-1):
                e = self.check_number(t[y][x+1] / t[y][x])
                if e == 0 or e > 1e9 or e is False:
                    return False
                seq.append(e)
            t.append(seq)
        return [KITAS(t, self.modulo_range, self.dimension)]

    def reduction(self, height):
        if len(self.eos) - height >= 2:
            t = KITAS(self.eos[height], self.modulo_range, self.dimension)
            return [t.basic()]
        return False

    def exponentiation(self, j):
        if len(self.eos[0]) < 2 or j == 0 or (0 in self.eos[0] and j <= 0):
            return False
        seq = []
        for i in self.eos[0]:
            e = (i)**(j)
            if type(e) == complex:
                return False
            e = self.check_number(e)
            if e is False:
                return False
            seq.append(e)
        seq = KITAS(seq, self.modulo_range, self.dimension)
        return [seq.basic()]

    def logarithm(self):
        base, ln, seq = self.eos[0], len(self.eos[0]), []
        if ln < 2 or 0 in base or 1 in base or -1 in base:
            return False
        for i in range(len(self.eos)-1):
            e = self.check_number(log(fabs(base[i+1]), fabs(base[i])))
            if e is False:
                return False
            seq.append(e)
        seq = KITAS(seq, self.modulo_range, self.dimension)
        return [seq.basic()]

    def double_operation(self, op1, op2):
        base, ln = self.eos[0], len(self.eos[0])
        seq, r = [], ln % 2
        cond1 = op1 == '/' and 0 in base[::2]
        cond2 = op2 == '/' and 0 in base[1::2]
        cond3 = op1 == '*' and r == 0 and 0 in base[:-2]
        cond4 = op2 == '*' and r != 0 and 0 in base[:-2]
        if ln < 2 or cond1 or cond2 or cond3 or cond4:
            return False
        for i in range(len(self.eos) - 1):
            op = op1 if i % 2 == 0 else op2
            e = self.do_op(op, base[i+1], base[i], False)
            if e is False:
                return False
            seq.append(e)
        edk_sol = KITAS(seq, self.modulo_range, self.dimension)
        return [edk_sol.basic()]

    def multi_level(self, dx, dy):
        if len(self.eos[0]) < 2:
            return False
        edks, p, x, y, dxy = [], 0, 0, 0, dx+dy
        for i in range(dxy):
            seq = []
            while len(self.eos[y]) > x and y < len(self.eos):
                seq.append(self.eos[y][x])
                if y+dy > len(self.eos) - 1:
                    break
                x, y = x+dx, y+dy
                p = p+1 if x == len(self.eos[y]) else p
            if len(seq) == 1 or (p == 0 and i == dxy-1 and dxy != 1):
                return False
            edk = KITAS(seq, self.modulo_range, self.dimension)
            edks.append(edk.basic())
            x, y = 0, i+1
        return edks

    def focusing(self, shift, divs):
        ln1, s = len(self.eos[0]), sum(divs)
        if shift >= s or s+shift > ln1 or ln1 < 4:
            return False
        main_list = self.eos[0][shift:]
        remainder, ln2 = self.eos[0][:shift], len(divs)
        sequences, edks = [[] for k in range(ln2)], []
        while main_list:
            for p in range(ln2):
                sequences[p] += main_list[:divs[p]]
                del main_list[:divs[p]]
                if not main_list:
                    break
        for q in range(ln2-1, -1, -1):
            if remainder:
                sequences[q] = remainder[-divs[q]:] + sequences[q]
                del remainder[-divs[q]:]
            if len(sequences[q]) <= 1:
                return False
            edk = KITAS(sequences[q], self.modulo_range, self.dimension)
            edks.insert(0, edk.basic())
        return edks

    def analogy(self, shift, nels):
        ln0 = len(self.eos[0][shift:])
        if ln0 < 4 or nels == 1 or nels+1 > ln0 or ln0 % nels == 0:
            return False
        seq0, sequences = deepcopy(self.eos[0][shift:]), []
        while seq0:
            seq = seq0[:nels] if len(seq0[:nels]) > 1 else [seq0[:nels]]
            sequences.append(KITAS(seq, self.modulo_range, self.dimension))
            del seq0[:nels]
        edks = list(map(lambda s: s.basic() if len(s.eos) > 1 else s, sequences))
        edkf = edks.pop()
        for i in range(1, len(edks[0].eos[0])):
            prz = list(map(lambda e: e.eos[i] == [0]*len(e.eos[i]), edks))
            nt = prz.count(True)
            if nt == len(prz):
                lnf = len(edkf.eos)
                if lnf == 1 or (i >= lnf and edkf.eos[-1] != [0]*len(edkf.eos[-1])):
                    h, i = len(edks[0].eos)-lnf, 1 if lnf == 1 else i
                    subedk = deepcopy(edks[0].eos[-h:])
                    if list(filter(lambda x: x.eos[-h:] != subedk, edks)):
                        return False
                    edkf.eos += subedk
                elif i < lnf and edkf.eos[i] == [0]*len(edkf.eos[i]) \
                        and edkf.eos[i-1] != [0]*len(edkf.eos[i-1]):
                    edkf.eos[i:] = deepcopy(edks[0].eos[i:])
                else:
                    return False
                if len(edkf.eos) != len(edks[0].eos):
                    return False
                ln2, ln1 = len(edkf.eos[i]), len(edkf.eos[i-1])-1
                ln3 = len(edkf.eos[0])
                SEQ = self.eos[0] + edkf.rev_basic(ln2-ln1, ln1, i)[ln3:]
                if SEQ == self.eos[0]:
                    return False
                return [KITAS([SEQ, 0], self.modulo_range, self.dimension)]
            elif nt != 0:
                return False
        return False

    def split_of_elements(self):
        check = list(map(lambda x: type(x) == float or x < 0, self.eos[0]))
        if len(self.eos[0]) <= 2 or True in check:
            return False
        s1, s2 = [], []
        for e in self.eos[0]:
            num_to_str = str(e)
            add_commas = ','.join(num_to_str)
            divide = add_commas.split(',')
            lend = [int(e) for e in divide]
            if lend != [lend[0]]*len(lend):
                return False
            s1.append(lend[0])
            s2.append(len(lend))
        if len(s1) < 2 or s1 == self.eos[0]:
            return False
        e1 = KITAS(s1, self.modulo_range, self.dimension)
        e2 = KITAS(s2, self.modulo_range, self.dimension)
        return [e1.basic(), e2.basic()]

    def dgee(self, x, y):
        if len(self.eos[0]) <= 2:
            return False
        s1, groups, a = [], [], None
        for e in self.eos[0]:
            if e != a:
                a = e
                s1.append(e)
                groups.append([])
            groups[-1].append(e)
        if len(s1) < 2:
            return False
        s2 = list(map(lambda y: len(y), groups))
        e1 = KITAS(s1[:-1] if x == 1 else s1, self.modulo_range, self.dimension)
        e2 = KITAS(s2[:-1] if y == 1 else s2, self.modulo_range, self.dimension)
        if len(e2.eos) < 2 or len(e1.eos) < 2 or s2 == [1]*len(s2):
            return False
        return [e1.basic(), e2.basic()]

    def dgde(self, x, y):
        seq, ln0 = self.eos[0][:], len(self.eos[0])
        if ln0 <= 2:
            return False
        s1 = [e for i, e in enumerate(seq) if e not in seq[:i]]
        ln1, u, v, c = len(s1), 0, 0, 0
        if self.eos[0][0:ln1] == s1:
            c = 1
        else:
            for j in range(ln1, 0, -1):
                if s1[:j] == seq[-len(s1[:j]):]:
                    u = len(s1[:j])
                    break
            if u == 0:
                return False
        s, groups = [self.eos[0].count(e) for e in s1], []
        while seq:
            g, ln2, ln3 = [], len(seq), u if c == 0 and v == 0 else ln1
            for i in range(ln3):
                if s[i] != 0:
                    g.append(s1[i])
                    s[i] = s[i]-1
            v, a = 1, (0, len(g)) if c == 1 else (ln2-len(g), ln2+1)
            if seq[a[0]:a[1]] != g:
                return False
            del seq[a[0]:a[1]]
            groups.append(g)
        groups = groups if c == 1 else groups[::-1]
        s2 = list(map(lambda z: len(z), groups))
        e1 = KITAS(s1[:-1] if x == 1 else s1, self.modulo_range, self.dimension)
        e2 = KITAS(s2[:-1] if y == 1 else s2, self.modulo_range, self.dimension)
        if len(e2.eos) < 2 or (len(e1.eos) < 2 and c == 0):
            return False
        return [e2.basic()] if c == 1 else [e1.basic(), e2.basic()]

    def repetition_symmetry(self, n):
        ln = len(self.eos[0])
        if ln <= 2:
            return False
        for i in range(2, ln):
            seq = deepcopy(self.eos[0])
            groups = [seq[k:k+i] for k in range(0, len(seq), i)]
            g0, gf, ln1 = groups[0], groups[-1], len(groups)
            lg0, lgf = len(g0), len(gf)
            if lg0 > lgf:
                for k in range(1, ln1):
                    if k == ln1-1 and g0[:lgf] == gf:
                        x = ceil((len(gf)+n)/len(g0))
                        s = [n for m in groups[:ln1-1] for n in m]+g0*x
                        return [KITAS([s, 0], self.modulo_range, self.dimension)]
                    if groups[k-1] != groups[k]:
                        break
        return False

    def specular_symmetry(self):
        ln = len(self.eos[0])
        if ln <= 2:
            return False
        for i in range(round(ln/2), ln):
            e0, e1 = self.eos[0][:i], self.eos[0][:i+1]
            e2, e3, e4 = self.eos[0][i:], e0[::-1], e1[::-1]
            if e0 == e2[::-1]:
                return False
            if e3[:len(e2)] == e2:
                if len(e0+e3) <= ln:
                    continue
                return [KITAS([e0+e3, 0], self.modulo_range, self.dimension)]
            if e4[:len(e2)] == e2 and i != ln-1:
                if len(e1+e4[1:]) <= ln:
                    continue
                return [KITAS([e1+e4[1:], 0], self.modulo_range, self.dimension)]
        return False

    def circular_symmetry(self, n):
        m, length, mask = 0, len(self.eos)+n, []
        if length < 3:
            return False
        for i in range(1, length+1):
            if i % 2 == 0:
                mask.append(m)
                m += 1
            else:
                mask.insert(0, m)
        ln = length if length % 2 != 0 else int(length/2)
        seq = self.eos[0] + [None]*n
        sol0 = self.check_csym(mask, ln, seq)
        if length % 2 == 0 and sol0 is False:
            m1 = list(range(1, ln))
            mask = [0] + m1 + [0] + m1[::-1]
            sol0 = self.check_csym(mask, ln, seq)
        return False if not sol0 else [KITAS([sol0, 0], self.modulo_range, self.dimension)]

    @staticmethod
    def check_csym(mask, length, seq):
        for j in range(length):
            ln = len(mask)
            link = list(zip(mask, seq))
            link0 = link[:]
            link.sort(key=lambda y: y[0])
            ln1 = ln if ln % 2 == 0 else ln-1
            sym = sum([1 for k in range(0, ln1, 2) if link[k][1] ==
                       link[k+1][1] or (link[k+1][1] is None and
                       link[k][1] is not None) or (link[k][1] is None
                       and link[k+1][1] is not None)])
            if sym == ln1 // 2:
                for k in range(len(link)):
                    if link[k][1] is None:
                        if link[k-1][0] != link[k][0] and k == len(link)-1:
                            continue
                        link0[link0.index(link[k])] = link[k-1] if link[k-1][0] == link[k][0] else link[k+1]
                sol = list(zip(*link0))
                if None not in sol[1]:
                    return list(sol[1])
            else:
                shift = mask.pop(-1)
                mask.insert(0, shift)
        return False

    def array_symmetry(self):
        length, dimension0, dimension1 = len(self.eos), self.dimension[0], self.dimension[1]
        if dimension0 <= 1 or dimension1 <= 1 or length < (dimension0-1)*dimension1+1 or length >= dimension0*dimension1:
            return False
        seq0, ln = self.eos[0][:], len(self.eos[0])
        seq0 += [None]*(dimension0*dimension1-ln)
        sr = self.rc_sym(seq0[:], dimension0, dimension1)
        if sr:
            return sr
        seq1 = []
        for i in range(dimension1):
            seq1 += seq0[i::dimension1]
        sc = self.rc_sym(seq1[:], dimension1, dimension0)
        if sc:
            sc, sct = sc[0].eos[0], []
            for i in range(dimension0):
                sct += sc[i::dimension0]
            return [KITAS([sct, 0], self.modulo_range, self.dimension)]
        if dimension0 != dimension1 or dimension0 <= 2:
            return False
        for j, k in zip([0, dimension0-1], [1, -1]):
            md = [i for i in range(j, len(seq0), dimension0+k)]
            seq2 = self.diagonal_sym(seq0[:], md[:], 0, -k)
            if seq2 is False:
                continue
            sd = self.diagonal_sym(seq2, md[:], -1, k)
            if sd is False:
                continue
            return KITAS([sd, 0], self.modulo_range, self.dimension)
        return False

    def rc_sym(self, seq, dimension0, dimension1):
        groups = [seq[p:p+dimension1] for p in range(0, len(seq), dimension1)]
        ln, sol = len(groups), []
        j = int(ln/2)
        i = j if ln % 2 != 0 else j-1
        for k in range(round(dimension0/2)):
            if None in groups[j]:
                d = groups[j].index(None)
                if groups[i][:d] != groups[j][:d]:
                    return False
                groups[j] = groups[i]
            else:
                if groups[i] != groups[j]:
                    return False
            i, j = i-1, j+1
        for m in groups:
            sol += m
        return [KITAS([sol, 0], self.modulo_range, self.dimension)]

    def diagonal_sym(self, seq, md, e1, e2):
        for j in range(self.dimension[0]-1):
            if j != 0:
                del md[e1]
                md = [k+e2 for k in md]
            md1 = [seq[m] for m in md]
            ln = len(md1)
            x, y = md1[:round(ln/2)], md1[ln//2:]
            if (None not in y and x != y[::-1]) or (None in y and None in x):
                return False
            if None in y:
                p, q = x[::-1], y.index(None)
                if p[:q] != y[:q]:
                    return False
                if len(x) == 1 or p[:q] == y[:q]:
                    md1[ln//2:] = p
                    for r in range(ln):
                        seq[md[r]] = md1[r]
        seq[-1] = seq[0] if seq[-1] is None else seq[-1]
        seq[0] = seq[-1] if seq[0] is None else seq[0]
        return seq

    def diagonal(self, a, b, c, d, e):
        dimension0, dimension1 = self.dimension[0], self.dimension[1]
        if dimension0 != dimension1-a or dimension0 <= b or dimension1 <= 2 or len(self.eos[0]) != dimension0*dimension1-c:
            return False
        array = self.eos[0]
        diagonal = [array[i] for i in range(d, len(array), e)]
        sol = KITAS(diagonal, self.modulo_range, self.dimension)
        return [sol.basic()]

    def reck(self, roc, ops):
        dimension0, dimension1 = self.dimension[0], self.dimension[1]
        ln, seq = len(ops), self.eos[0][:]
        if dimension0*dimension1-1 != len(seq) or (roc == 'c' and (dimension0 < 3 or dimension1 < 2)) or (roc == 'r' and (dimension0 < 2 or dimension1 < 3)):
            return False
        groups = [seq[dimension1*j:dimension1*(j+1)] for j in range(dimension0)] if roc == 'r' else [seq[j::dimension1] for j in range(dimension1)]
        o, dimensionx = ('/','*','-','+'), dimension1-2 if roc == 'r' else dimension0-2
        combs, sol = tuple(product(ops, repeat=dimensionx)), None
        for c in combs:
            pos_sol = self.reck2(deepcopy(groups), list(c), c, o, dimensionx)
            if pos_sol:
                sol = pos_sol
                break
        return False if sol is None else [KITAS([self.eos[0]+[sol], 0], self.modulo_range, self.dimension)]

    def reck2(self, groups, ops, c, o, dimension):
        for g in groups:
            e = g.pop() if g != groups[-1] else None
            m = [0]*(len(g)+dimension)
            m[::2], m[1::2] = g, ops
            for j in range(4):
                for k in range(c.count(o[j])):
                    i = m.index(o[j])
                    x = self.do_op(o[j], m[i-1], m[i+1], False)
                    if x is False:
                        return False
                    m[i-1:i+2] = [x]
            if m[0] != e and e is not None:
                return False
            elif e is None:
                return m[0]

    def transposed(self):
        dimension0, dimension1 = self.dimension[0], self.dimension[1]
        if dimension0 <= 1 or dimension1 <= 1 or len(self.eos[0]) != dimension0*dimension1-1:
            return False
        transposed = []
        for i in range(dimension1):
            transposed += self.eos[0][i::dimension1]
        sol = KITAS(transposed, self.modulo_range, (dimension1, dimension0))
        return [sol.basic()]
    
    def rev_basic(self, n, m, h):
        for j in range(h, 0, -1):
            for k in range(n):
                if k+m < len(self.eos[j]):
                    e = self.check_number(self.eos[j][k+m]+self.eos[j-1][k+m])
                    if e is False:
                        return False
                    self.eos[j-1].append(e)
            m += 1
        return self.eos[0]

    def rev_divisions(self, n, m, h):
        for j in range(h, 0, -1):
            for k in range(n):
                if k+m < len(self.eos[j]):
                    e = self.check_number(self.eos[j][k+m]*self.eos[j-1][k+m])
                    if e is False:
                        return False
                    self.eos[j-1].append(e)
            m += 1
        return self.eos[0]

    def rev_exponentiation(self, j):
        if j == -1 and 0 in self.eos[0]:
            return False
        row = []
        for x in self.eos[0]:
            try:
                e = (x)**(1/j)
                if type(e) == complex:
                    return False
                e = self.check_number(e)
                if e is False:
                    return False
                row.append(e)
            except:
                row.append(inf)
        self.eos[0] = row
        return self.eos[0]

    def rev_logarithm(self, n, m):
        for k in range(n):
            if k+m < len(self.eos[1]):
                if self.eos[0][k+m] == 0 and self.eos[1][k+m] == 0:
                    return False
                try:
                    e = (self.eos[0][k+m])**(self.eos[1][k+m])
                    if type(e) == complex:
                        return False
                    e = self.check_number(e)
                    if e is False:
                        return False
                    self.eos[0].append(e)
                except:
                    self.eos[0].append(inf)
        return self.eos[0]

    def rev_double_operation(self, op1, op2, n, m):
        ln_2 = len(self.eos[0]) % 2
        for k in range(n):
            if k+m < len(self.eos[1]):
                op = op1 if (ln_2 == 0 and k % 2 != 0) or (ln_2 != 0 and k % 2 == 0) else op2
                e = self.do_op(op, self.eos[1][k+m], self.eos[0][k+m], True)
                if e is False:
                    return False
                self.eos[0].append(e)
        return self.eos[0]

    def rev_multi_level_p(self, div):
        base_copy = self.eos[0][:]
        for j in range(len(self.eos) - 1, 0, -1):
            for k in range(len(self.eos[j])):
                if k >= len(self.eos[j-1]):
                    continue
                if self.eos[j][k] is not None and self.eos[j-1][k] is not None:
                    e0 = self.eos[j][k]*self.eos[j-1][k] if div else self.eos[j][k]+self.eos[j-1][k]
                    e = self.check_number(e0)
                    if e is False:
                        return False
                    if k < len(self.eos[j-1]) - 1:
                        self.eos[j-1][k+1] = e
                    else:
                        self.eos[j-1].append(e)
        seq = self.eos[0]
        if None in seq:
            seq = self.eos[0][:self.eos[0].index(None)]
        return seq if seq != base_copy else False

    def rev_multi_level_i(self, subedk, dx, dy, b):
        ln, vx, vy = len(subedk.eos[0]), [], []
        for i in range(ln):
            vx.append(dx*i)
            vy.append(dy*i+b)
        if dx == 0 or dy == 0:
            x, y = ln, ln
        else:
            x, y = dx*ln, dy*ln+b
        for j in range(y):
            if j > len(self.eos)-1:
                self.eos.append([None]*x)
            if x > len(self.eos[j]):
                self.eos[j] += [None]*(x-len(self.eos[j]))
        for k in range(ln):
            p, q = vx[k], vy[k]
            if self.eos[q][p] is None:
                self.eos[q][p] = subedk.eos[0][k]
        return self.eos

    def rev_focusing_i(self, subedks, shift, divs, n, ln):
        seq = self.eos[0]+[None]*n
        for i in range(ln):
            l = subedks[i] if type(subedks[i]) == list else subedks[i].eos[0]
            a = shift + sum(divs[:i])
            s, o = sum(divs), round(len(l)/divs[i])+1
            race = [a+j+s*k for k in range(o) for j in range(divs[i]) if a+j+s*k <= len(seq)-1]
            for p, q in zip(race, range(len(race))):
                if q < len(l) and seq[p] is None:
                    seq[p] = l[q]
        return seq if None not in seq else seq[:seq.index(None)]

    def rev_split_of_elements_i(self, subedks):
        s1, s2, r = subedks[0].eos[0], subedks[1].eos[0], []
        for x, y in zip(s1, s2):
            if y == 0:
                return False
            if (type(x) == float and x-int(x) != 0.0) or x < 0 or \
                    (type(y)==float and y-int(y) != 0.0) or y < 0 or y > 10**5:
                return False
            r.append(int(str(int(x))*int(y)))
        return r if self.eos[0] == r[:len(self.eos[0])] else False

    def rev_dgee_i(self, subedks):
        s1, s2, r = subedks[0].eos[0], subedks[1].eos[0], []
        try:
            for x, y in zip(s1, s2):
                if y < 0 or (type(y) != int and y-int(y) != 0.0):
                    return False
                if y > 200:
                    break
                r += [x]*int(y)
                if len(r) > 200:
                    break
            return r if self.eos[0] == r[:len(self.eos[0])] else False
        except:
            return False

    def rev_dgde_i(self, subedks):
        if len(subedks) == 1:
            s1 = []
            for e in self.eos[0]:
                if e not in s1:
                    s1.append(e)
            s2 = subedks[0].eos[0]
        else:
            s1 = subedks[0].eos[0]
            s2 = subedks[1].eos[0]
        seq, r1, r2 = self.eos[0][:], [], []
        for i in s2:
            if i < 0 or (type(i) != int and i-int(i) != 0.0):
                return False
            i = int(i)
            r1, r2 = r1+s1[:i], r2+s1[-i:]
            del seq[:i]
        if self.eos[0] == r1[:len(self.eos[0])]:
            return r1
        elif self.eos[0] == r2[:len(self.eos[0])]:
            return r2
        else:
            return False