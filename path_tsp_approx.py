import pulp
import itertools
from graph import *
import copy
from typing import *


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s) + 1))


class PathTSPProblem:
    """
    Represents an s-t path TSP problem.
    """

    def __init__(self, graph, s, t):
        """
        graph must be metricgraph
        """
        self.s = s
        self.t = t
        self.graph = graph


def delta(nodeset: Tuple[Node,...], edges: List[Edge], emptyindex: int = None) -> List[int]:
    l = []
    for i, edge in enumerate(edges):
        if len({edge.n1, edge.n2}.intersection(set(nodeset))) == 1:
            l.append(i)
    if len(nodeset) == 0 and emptyindex is not None:
        return [emptyindex]
    return l


def held_karp_lp(pathtsp: PathTSPProblem, integer: int = False, b_cuts: Set[List[int]] = None):
    """
    returns lp problem representing held-karp
    """

    prob: pulp.LpProblem = pulp.LpProblem('held-karp', sense=pulp.LpMinimize)

    # we have one variable per edge in the complete graph
    lp_x: List[pulp.LpVariable] = []
    edges = []  # for keeping track of index
    edges_to_index: Dict[Edge, int] = {}  # convention: smallest then biggest label
    for e in pathtsp.graph.edges():
        edges_to_index[e] = len(lp_x)
        edges.append(e)
        if not integer:
            lp_x.append(pulp.LpVariable(str(e), lowBound=0))
        else:
            lp_x.append(pulp.LpVariable(str(e), lowBound=0, cat='Integer'))

    # objective (eh.... how to create empty thing here)
    objective: pulp.pulp.LpAffineExpression = edges[0].w * lp_x[0]
    for i, edge in enumerate(edges):
        if i == 0:
            continue
        objective += edges[i].w * lp_x[i]
    prob += objective

    def delta_part(nodeset: Tuple[Node,...], edges: List[Edge]) -> pulp.pulp.LpAffineExpression:
        d_i = delta(nodeset, edges)
        f = lp_x[d_i[0]]
        for d in d_i[1:]:
            f += lp_x[d]
        return f

    # now add all constraints
    # we do it in the order that held-karp is presented on page 1540 of 1.5-approx paper
    # iterate over all subsets of nodes
    for nodeset in powerset(pathtsp.graph.nodes()):
        # figure out which type of set this is
        stcount = 0
        if pathtsp.s in nodeset:
            stcount += 1
        if pathtsp.t in nodeset:
            stcount += 1
        if len(nodeset) > 0 and len(nodeset) < len(pathtsp.graph.nodes()) and stcount in {0, 2}:
            # case 1!!!
            # generate edges
            # do it naively becuase meh
            prob += delta_part(nodeset, edges) >= 2
        elif stcount == 1:
            # case 2!!
            prob += delta_part(nodeset, edges) >= 1

        if len(nodeset) == 1 and stcount == 0:
            # case 3!
            prob += delta_part(nodeset, edges) == 2
        elif len(nodeset) == 1 and stcount == 1:
            # case 4!
            prob += delta_part(nodeset, edges) == 1

    # if we have any b_cuts, also add those constraints
    if b_cuts is not None:
        for nodeset in b_cuts:
            prob += delta_part(nodeset, edges) >= 3

    return prob, {e: x for e, x in zip(edges, lp_x)}


def find_optimal_tsp(pathtsp):
    """
    return the optimal path tsp solution, as a list of edges
    """
    x_star = find_x_star(pathtsp, integer=True)

    # just take all edges used
    opt = []
    for e in x_star:
        for num in range(round(x_star[e])):
            opt.append(e)
    return opt


# this function finds x^*
def find_x_star(pathtsp:PathTSPProblem, integer:int=False, b_cuts:Set[List[int]]=None)->Dict[Edge,float]:
    """
    pathtsp: tsp problem
    finds x^*; returns as dictionary edge -> value
    """

    hk, edges = held_karp_lp(pathtsp, integer=integer, b_cuts=b_cuts)

    status = hk.solve()

    # print(hk)
    # print(pulp.LpStatus[status])
    # print(pulp.value(hk.objective))
    # print(hk.objective)
    # print([pulp.value(x) for x in hk.variables()])
    # print([str(e) for e in edges])

    x_star = {e: pulp.value(edges[e]) for e in edges}

    return x_star


def find_B(pathtsp:PathTSPProblem, x)->Set[List[int]]:
    """
    returns a list of cuts: all cuts in B(x)
    """
    B = set()
    # iterate over all subsets
    edges = pathtsp.graph.edges()
    for nodeset in powerset(pathtsp.graph.nodes()):
        # only cuts with s in it and t not in it
        if pathtsp.s not in nodeset:
            continue
        if pathtsp.t in nodeset:
            continue

        # only cuts such that x(delta(nodeset)) < 3
        this_cut_value = 0
        for i in delta(nodeset, edges):
            this_cut_value += x[edges[i]]

        if this_cut_value < 3:
            B.add(nodeset)
    return B


INFTY = 100000000000


def get_y(pathtsp, B_family):
    # create a new family
    new_b = list(copy.deepcopy(B_family))
    new_b.append(tuple())  # the empty set
    final_b = tuple(pathtsp.graph.nodes())
    value, y_star = get_y_recursive(pathtsp, new_b, final_b, pathtsp.t)
    return y_star


def get_y_recursive(pathtsp, B_family, B, t, dp_memo=None):
    """
    finds y recursively, by trying all possible values for where
    the previous 1-cut could be
    """
    if dp_memo is None:
        dp_memo = {}
    if (B, t) in dp_memo:
        return dp_memo[(B, t)]
    print(f"Now considering cut {B} with end node {t} for the first time")
    if len(B) == 0:
        # base case!!!
        return 0, {}
    edges = pathtsp.graph.edges() + [Edge(Node('dummy'), pathtsp.s, 0)]
    minvalhere = INFTY
    minsolhere = None
    for B_prime in B_family:
        # check if B_prime is subset of B
        if not set(B_prime) < set(B):
            continue
        if t in B_prime:
            continue
        for ei in delta(B_prime, edges, emptyindex=len(edges) - 1):
            e_prime = edges[ei]
            # print(f'B: {B} B_prime: {B_prime} t: {t} e_prime: {e_prime}')
            # make sure we have a valid chain
            intersectgraph = set(B).difference(set(B_prime))
            if len({e_prime.n1, e_prime.n2}.intersection(intersectgraph)) != 1:
                continue
            e_1 = e_prime.n1
            e_2 = e_prime.n2
            if e_1 in intersectgraph:
                e_1 = e_prime.n2
                e_2 = e_prime.n1
            # now call dp on this
            this_val, this_sol = get_y_recursive(pathtsp, B_family, B_prime, e_1, dp_memo=dp_memo)
            # now solve the LP
            # we need to be careful here because we might have that the start node
            # is equal to the end node. that's fine tho.
            if len(intersectgraph) == 1:
                if e_2 != t:
                    raise RuntimeError('something very very wrong')
                # this means we can just use thisAvalue, plus the edge
                if this_val + e_prime.w < minvalhere:
                    minvalhere = this_val + e_prime.w
                    minsolhere = copy.deepcopy(this_sol)
                    print(minsolhere)
                    if e_prime in minsolhere:
                        raise RuntimeError('should never happen')
                    minsolhere[e_prime] = 1
            else:
                # now we need to solve the LP!!
                # exciting
                # so what do we do here? create a new metric graph
                # so we want to create a graph defined on a subset
                induced_graph = pathtsp.graph.subset_graph(intersectgraph)
                # print(intersectgraph)
                # print(e)
                # print(B_prime)
                # print(e_prime)
                sub_pathtsp = PathTSPProblem(induced_graph, e_2, t)
                # print(induced_graph)
                # print(pathtsp.graph)
                relevant_b_cuts = [set(b).intersection(intersectgraph) for b in B_family if
                                   set(b) < set(B) and set(B_prime) < set(b)]
                y_i = find_x_star(sub_pathtsp, b_cuts=relevant_b_cuts)

                val_y_i = sum([k.w * y_i[k] for k in y_i])

                if this_val + val_y_i + e_prime.w < minvalhere:
                    minvalhere = this_val + val_y_i + e_prime.w
                    minsolhere = copy.deepcopy(this_sol)
                    print(minsolhere)
                    for k in y_i:
                        if k in minsolhere:
                            raise RuntimeError('should never happen')
                        minsolhere[k] = y_i[k]
                    if e_prime in minsolhere:
                        raise RuntimeError('should never happen')
                    minsolhere[e_prime] = 1

    dp_memo[(B, t)] = (minvalhere, copy.deepcopy(minsolhere))
    return dp_memo[(B, t)]


def integrality_gap_ptsp_graph(n):
    # look at page 2 of shmoys
    v = [Node('s'), Node('t')] + [Node(f't{i}') for i in range(n)] + [Node(f'b{i}') for i in range(n)]
    # t is top row b is bottom row
    e = [Edge(Node('s'), Node('t0'), 1), Edge(Node('s'), Node('b0'), 1), Edge(Node('t0'), Node('b0'), 1)] + \
        [Edge(Node('t'), Node(f't{n - 1}'), 1), Edge(Node('t'), Node(f'b{n - 1}'), 1),
         Edge(Node(f't{n - 1}'), Node(f'b{n - 1}'), 1)] + \
        [Edge(Node(f't{i}'), Node(f't{i + 1}'), 1) for i in range(n - 1)] + \
        [Edge(Node(f'b{i}'), Node(f'b{i + 1}'), 1) for i in range(n - 1)]
    s = SimpleGraph(v, e)

    s.to_dot('ssss.dot')

    m = s.to_metric()
    return m, e


def integrality_gap_ptsp_graph_perturbed(n):
    # look at page 2 of shmoys
    v = [Node('s'), Node('t')] + [Node(f't{i}') for i in range(n)] + [Node(f'b{i}') for i in range(n)]
    # t is top row b is bottom row
    e = [Edge(Node('s'), Node('t0'), 1), Edge(Node('s'), Node('b0'), 0.3), Edge(Node('t0'), Node('b0'), 1)] + \
        [Edge(Node('t'), Node(f't{n - 1}'), 1 - 0.3), Edge(Node('t'), Node(f'b{n - 1}'), 1),
         Edge(Node(f't{n - 1}'), Node(f'b{n - 1}'), 1)] + \
        [Edge(Node(f't{i}'), Node(f't{i + 1}'), 1 + 0.2) for i in range(n - 1)] + \
        [Edge(Node(f'b{i}'), Node(f'b{i + 1}'), 1) for i in range(n - 1)]
    s = SimpleGraph(v, e)

    s.to_dot('ssss.dot')

    m = s.to_metric()
    return m, e


def y_not_0_graph(n):
    # look at page 2 of shmoys
    v = [Node('s'), Node('t')] + [Node(f't{i}') for i in range(n)] + [Node(f'b{i}') for i in range(n)] + \
        [Node(f'tt{i}') for i in range(n)] + [Node(f'bb{i}') for i in range(n)] + \
        []
    # [Node('s1'), Node('s2')]
    # t is top row b is bottom row
    e = [Edge(Node('s'), Node('t0'), 1), Edge(Node('s'), Node('b0'), 1), Edge(Node('t0'), Node('b0'), 1)] + \
        [Edge(Node('t'), Node(f't{n - 1}'), 1), Edge(Node('t'), Node(f'b{n - 1}'), 1),
         Edge(Node(f't{n - 1}'), Node(f'b{n - 1}'), 1)] + \
        [Edge(Node(f't{i}'), Node(f't{i + 1}'), 1) for i in range(n - 1)] + \
        [Edge(Node(f'b{i}'), Node(f'b{i + 1}'), 1) for i in range(n - 1)] + \
        [Edge(Node(f'tt{i}'), Node(f'tt{i + 1}'), 1) for i in range(n - 1)] + \
        [Edge(Node(f'bb{i}'), Node(f'bb{i + 1}'), 1) for i in range(n - 1)] + \
        [Edge(Node('s'), Node('tt0'), 1), Edge(Node('s'), Node('bb0'), 1), Edge(Node('tt0'), Node('bb0'), 1)] + \
        [Edge(Node('t'), Node(f'tt{n - 1}'), 1), Edge(Node('t'), Node(f'bb{n - 1}'), 1),
         Edge(Node(f'tt{n - 1}'), Node(f'bb{n - 1}'), 1)] + \
        []
    # [Edge(Node('s1'), Node('b0'), 1), Edge(Node('s1'),Node('b3'), 1)] + \
    # [Edge(Node('s2'), Node('t0'), 1), Edge(Node('s2'),Node('t3'), 1)] + \
    # [Edge(Node('s1'),Node('b1'),1), Edge(Node('s1'),Node('b2'), 1)]
    s = SimpleGraph(v, e)

    s.to_dot('ssss.dot')

    m = s.to_metric()
    return m, e


def k_branch_graph(k, n):
    """
    create a k-branch graph
    """
    v = [Node('s'), Node('t')]

    for c in range(k):
        for i in range(n):
            v.append(Node(f'{c}x{i}'))

    e = []

    # connect s with 0 nodes in full clique
    left_nodes = [Node('s')] + [Node(f'{c}x0') for c in range(k)]
    for n1, n2 in itertools.combinations(left_nodes, 2):
        e.append(Edge(n1, n2, 1))  # use weight 1 for now

    # connect t with n-1 nodes in full clique
    right_nodes = [Node('t')] + [Node(f'{c}x{n - 1}') for c in range(k)]
    for n1, n2 in itertools.combinations(right_nodes, 2):
        e.append(Edge(n1, n2, 1))  # use weight 1 for now

    # now connect each branch within itself
    for c in range(k):
        for i in range(n - 1):
            e.append(Edge(Node(f'{c}x{i}'), Node(f'{c}x{i + 1}'), 1))

    s = SimpleGraph(v, e)

    s.to_dot('ssss.dot')

    m = s.to_metric()
    return m, e


def final_try_graph(n):
    """
    create graph where hopefully choosing the spanning
    tree from the optimum LP solution is not very good
    create 3n+1 nodes
    """
    v = [Node(str(i)) for i in range(3 * n + 1)]

    eps = 0.1

    e = []
    # straight edges
    e += [Edge(Node(f'{i}'), Node(f'{i + 1}'), 1) for i in range(3 * n)]

    # jump-ahead 3 steps
    w = [1] + [1 - eps for i in range(n - 2)] + [1]
    e += [Edge(Node(f'{3 * i}'), Node(f'{3 * i + 3}'), w[i]) for i in range(n)]

    # jump ahead 2 steps
    e += [Edge(Node(f'{3 * i + 2}'), Node(f'{3 * i + 4}'), 1) for i in range(n - 1)]

    # final jump ahead for some reason
    e += [Edge(Node('1'), Node(f'{3 * n - 1}'), 1)]

    # i think this is it

    s = SimpleGraph(v, e)
    s.to_dot('ssss.dot')

    m = s.to_metric()
    return m, e
