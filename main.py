from timeit import default_timer
from math import exp
from kitas import KITAS
from search import SeqSearchAlgorithm,SeqState
from predict import SeqPredictor
class KitBit:
    def __init__(self,structure,kl,mni,depth,search_algorithm='BFS',n=3,min_zeros=1,epsilon=exp(-18),all_solutions=False):
        self.structure = structure
        self.kl = kl
        self.mni = mni
        self.depth = depth
        self.search_alg = search_algorithm
        self.n = n
        self.min_zeros = min_zeros
        self.epsilon = epsilon
        self.all_solutions = all_solutions

    def solve_seq(self):
        if len(self.structure) == 1:
            return (False, 0, self.kl)
        init_time = default_timer()
        sol = self.numeric_solver(self.structure[:], self.kl, self.depth, 10)
        end_time = default_timer()
        return [sol, end_time - init_time]

    def numeric_solver(self, seq, kl, depth, module):
        edk0 = KITAS(seq[:], module, None)
        edk0 = edk0.basic()
        if edk0.is_goal('BASIC', self.epsilon, self.min_zeros, 0):
            edk0.eos[-1] += [0]*self.n
            return edk0.rev_basic(self.n, 1, len(edk0.eos)-1), 0, ['BASIC']
        if len(seq) == 2:
            edk0.eos[-1] += edk0.eos[-1]*self.n
            return edk0.rev_basic(self.n, 1, len(edk0.eos)-1), 0, ['BASIC']
        st0 = SeqState(None, [edk0], [], [])
        alg = SeqSearchAlgorithm(st0, kl, self.mni, depth, self.n, self.min_zeros, self.epsilon)
        if self.search_alg == 'BFS':
            lst = alg.bfs()
        elif self.search_alg == 'DFS':
            lst = alg.dfs([st0])
        else:
            lst = alg.branch()
        if lst is None or lst is False or len(lst[0]) == 1:
            return False, 0, kl
        acts = [st.action for st in lst[0] if st.action is not None]
        sol = SeqPredictor(lst[0][::-1], self.n)
        ln, sol = len(seq), sol.predictor()
        return (False,0,kl) if not sol or len(sol)<=len(seq) else (seq[:]+sol[ln:ln+self.n],lst[1],acts)