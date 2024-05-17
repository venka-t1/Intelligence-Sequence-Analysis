from kitas import KITAS
from copy import deepcopy
class SeqSearchAlgorithm:
    def __init__(self, init_state, kl, mni, depth, n, min_zeros, epsilon):
        self.init_state = init_state
        self.kl = kl
        self.mni = mni
        self.depth = depth
        self.n = n
        self.min_zeros = min_zeros
        self.epsilon = epsilon
        self.road = [init_state]
        self.count = 0
        self.objective = 'None'

    def bfs(self):
        i, k, ln, self.road = 0, -1, len(self.kl)-1, [self.road]
        n_iter = [len(self.kl)**p for p in range(1, self.depth+1)]
        for j in range(1, sum(n_iter)+1):
            k = k+1 if k < ln else 0
            st = False if not self.road[i] else \
                self.road[i][-1].new_state(self.n, self.road[i][-1].action,
                    self.kl[k], self.epsilon, self.min_zeros)
            if st is False:
                i = i+1 if k == ln else i
                self.road.append(False)
                continue
            self.road.append(self.road[i]+[st])
            i = i+1 if k == ln else i
            if len(st.sols) == len(st.eos):
                return self.road[-1], j
            if j >= self.mni:
                return False
        return False

    def dfs(self, road=[], level=0):
        for i in range(len(self.kl)):
            self.count += 1
            if len(self.road) > 1 or not self.objective:
                break
            st = road[-1].new_state(self.n, road[-1].action, self.kl[i], self.epsilon, self.min_zeros)
            if st is False:
                continue
            road.append(st)
            if len(st.sols) == len(st.eos):
                self.road = deepcopy(road)
                return road, self.count
            if self.count >= self.mni:
                self.objective = False
                return False
            if level < self.depth - 1:
                self.dfs(road, level+1)
            road.pop()
        if self.road:
            return self.road, self.count

    def kb_seq_gb(self):
        branches, self.road = [], [self.road]
        for i in range(self.depth):
            branches1 = []
            for lst in self.road:
                sts = list(map(lambda kita: lst[-1].new_state(self.n, lst[-1].action, kita, self.epsilon, self.min_zeros), self.kl))
                self.count += len(sts)
                sts = list(filter(lambda x: x is not False, sts))
                sts_sol = list(filter(lambda x: len(x.sols) == len(x.eos), sts))
                sts_not_sol = list(filter(lambda x: len(x.sols) != len(x.eos), sts))
                branches1 += list(map(lambda st: deepcopy(lst)+[st], sts_not_sol))
                branches += list(map(lambda st: deepcopy(lst)+[st], sts_sol))
            self.road = branches1
            if not branches1:
                break
        return (branches, self.count) if branches else False

    def branch(self):
        ln, road = len(self.kl), self.road
        for j in range(ln):
            st = False if not self.road[j] else road[-1].new_state(self.n, road[-1].action, self.kl[j], self.epsilon, self.min_zeros)
            if st is False:
                return False
            road.append(st)
            if len(st.sols) == len(st.eos) and j == ln-1:
                return road, j
        return False

class SeqState:
    def __init__(self, action, eos, sols, parents):
        self.action = action
        self.eos = eos
        self.sols = sols
        self.parents = parents

    def new_state(self, n, pva, act, epsilon, min_zeros):
        ln, edks, parents = len(self.eos), [], []
        for i in range(ln):
            edk = KITAS(self.eos[i].eos, self.eos[0].modulo_range)
            if (ln > 1 and i not in self.sols) or ln == 1:
                try:
                    new_edk = self.do_action(edk, n, pva, act)
                except:
                    return False
                if new_edk is False or False in new_edk:
                    return False
                edks += new_edk
                parents += [i for q in range(len(new_edk))]
        m = 0 if act != 'DIV' else 1
        sols = [e for e in range(len(edks)) if edks[e].is_goal(act, epsilon, min_zeros, m)]
        return SeqState(act, edks, sols, parents)

    @staticmethod
    def do_action(edk, n, pva, act):
        if act == 'DIV':
            return edk.divisions()
        if act[:3] == 'RED':
            a = act[act.index('(')+1:act.index(')')]
            return edk.reduction(int(a))
        elif act[:3] == 'EXP':
            a = act[act.index('(')+1:act.index(')')]
            if '/' in a:
                e0 = a.split('/')
                e = int(e0[0]) / int(e0[1])
            elif '.' in a:
                e = float(a)
            else:
                e = int(a)
            return edk.exponentiation(e)
        elif act == 'LOG':
            return edk.logarithm()
        elif act[:3] == 'DOP':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.double_operation(a[0], a[1])
        elif act[:2] == 'ML':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.multi_level(int(a[0]), int(a[1]))
        elif act[:3] == 'FOC':
            a = act[act.index('(')+1:act.index(')')].split(',')
            divs = [int(i) for i in a[1:]]
            return edk.focusing(int(a[0]), divs)
        elif act[:3] == 'ANA':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.analogy(int(a[0]), int(a[1]))
        elif act[:3] == 'SOE':
            return edk.split_of_elements()
        elif act[:4] == 'DGEE':
            a = act[act.index('(') + 1:act.index(')')].split(',')
            return edk.dgee(int(a[0]), int(a[1]))
        elif act[:4] == 'DGDE':
            a = act[act.index('(') + 1:act.index(')')].split(',')
            return edk.dgde(int(a[0]), int(a[1]))
        elif act == 'RSYM':
            return edk.repetition_symmetry(n)
        elif act == 'SSYM':
            return edk.specular_symmetry()
        elif act == 'CSYM':
            return edk.circular_symmetry(n)
        elif pva is not None and ('DIAG' in pva or pva == 'LOG' or
                pva[:3] == 'RED' or pva[:2] == 'ML' or pva[:3] == 'DOP'
                or pva[:3] == 'FOC' or pva[:3] == 'SOE' or 'DG' in pva):
            return False
        elif act == 'ASYM':
            return edk.array_symmetry()
        elif act == 'LDIAG':
            return edk.diagonal(0, 2, 1, 0, edk.dimension1+1)
        elif act == 'RDIAG':
            return edk.diagonal(1, 1, 0, edk.dimension1-1, edk.dimension1-1)
        elif act[:4] == 'RECK':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.reck(a[0], a[1])
        else:
            return edk.transposed()