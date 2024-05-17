class SeqPredictor:
    def __init__(self, states_list, n):
        self.states_list = states_list
        self.n = n

    def predictor(self):
        for i in range(len(self.states_list)-1):
            st0, st1 = self.states_list[i-1], self.states_list[i]
            act0, act1, p = st0.action, st1.action, st1.parents
            for j in range(len(st1.eos)):
                edk = st1.eos[j]
                if j in st1.sols:
                    if act1 == 'DIV':
                        edk.eos[-1] += [1]*self.n
                        seq = edk.rev_divisions(self.n, 1, len(edk.eos)-1)
                    elif act1[:3] == 'ANA' or act1 == 'RECK' or 'SYM' in act1:
                        seq = st1.eos[j].eos[0]
                    else:
                        edk.eos[-1] += [0]*self.n
                        seq = edk.rev_basic(self.n, 1, len(edk.eos)-1)
                elif act0[:3] == 'RED' or act0 == 'LOG' or act0[:2] == 'ML' or act0[:3] == 'DOP':
                    seq = self.prediction(act0, act1, edk)
                else:
                    seq = st1.eos[j].eos[0]
                if act1[:3] == 'EXP':
                    if not seq:
                        return False
                    a = act1[act1.index('(')+1:act1.index(')')]
                    if '/' in a:
                        e0 = a.split('/')
                        e = int(e0[0]) / int(e0[1])
                    elif '.' in a:
                        e = float(a)
                    else:
                        e = int(a)
                    seq = edk.rev_exponentiation(e)
                if not seq:
                    return False
                new_edk = self.insertion(st1.eos, p, st1.action, self.states_list[i+1].eos[p[j]], seq, j)
                if new_edk is False:
                    return False
                self.states_list[i+1].eos[p[j]].eos = new_edk
        ef, af = self.states_list[-1].eos[0], self.states_list[-2].action
        return self.prediction(af, af, ef) if af[:3] == 'RED' or \
            af == 'LOG' or af[:2] == 'ML' or af[:3] == 'DOP' else ef.eos[0]

    def prediction(self, act0, act1, edk):
        if act0 == 'LOG':
            return edk.rev_logarithm(self.n, len(edk.eos)-1)
        elif act0[:3] == 'DOP':
            return edk.rev_double_operation(
                act0[4], act0[6], self.n, len(edk.eos)-1)
        elif act0[:3] == 'RED':
            h = int(act0[4])
            if act1[:3] == 'DIV':
                return edk.rev_divisions(self.n, len(edk.eos)-h, h)
            return edk.rev_basic(self.n, len(edk.eos)-h, h)
        return edk.rev_multi_level_p(act1[:3] == 'DIV')

    def insertion(self, edks1, p, act1, edk2, seq, j):
        if act1[:2] == 'ML':
            a, parents1 = [], p[:]
            while parents1:
                b = p.count(p[0])
                a = a + list(range(b))
                del parents1[:b]
            act1 = act1[act1.index('(')+1:act1.index(')')].split(',')
            b, dx, dy = a[j], int(act1[0]), int(act1[1])
            return edk2.rev_multi_level_i(edks1[j], dx, dy, b)
        elif act1[:3] == 'FOC':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                act1 = act1[act1.index('(')+1:act1.index(')')].split(',')
                shift, divs = int(act1[0]), [int(x) for x in act1[1:]]
                edk2.eos[0] = edk2.rev_focusing_i(edks1, shift, divs, self.n, p.count(p[j]))
        elif act1[:3] == 'SOE':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                edk2.eos[0] = edk2.rev_split_of_elements_i(edks1)
                return False if not edk2.eos[0] else edk2.eos
        elif act1[:4] == 'DGEE':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                edk2.eos[0] = edk2.rev_dgee_i(edks1)
                return False if not edk2.eos[0] else edk2.eos
        elif act1[:4] == 'DGDE':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                edk2.eos[0] = edk2.rev_dgde_i(edks1)
                return False if not edk2.eos[0] else edk2.eos
        elif act1[:3] == 'RED':
            edk2.eos[int(act1[4])] = seq
        elif act1[:3] == 'DOP' or act1 == 'LOG':
            edk2.eos[1] = seq
        elif 'DIAG' in act1 or act1 == 'TRA':
            edk2.eos[0] = seq+[seq[edk2.dimension[0]-1]]
        else:
            edk2.eos[0] = seq
        return edk2.eos