import pulp
import itertools
from graph import *
import copy

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

        

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


def delta(nodeset, edges, emptyindex=None):
    l = []
    for i, edge in enumerate(edges):
        if len({edge.n1,edge.n2}.intersection(set(nodeset))) == 1:
            l.append(i)
    if len(nodeset) == 0 and emptyindex is not None:
        return [emptyindex]
    return l


def held_karp_lp(pathtsp, integer=False, b_cuts=None):
    """
    returns lp problem representing held-karp
    """

    prob = pulp.LpProblem('held-karp', sense=pulp.LpMinimize)

    # we have one variable per edge in the complete graph
    lp_x = []
    edges = [] # for keeping track of index
    edges_to_index = {} # convention: smallest then biggest label
    for e in pathtsp.graph.edges():
        edges_to_index[e] = len(lp_x)
        edges.append(e)
        if not integer:
            lp_x.append(pulp.LpVariable(str(e), lowBound=0))
        else:
            lp_x.append(pulp.LpVariable(str(e), lowBound=0, cat = 'Integer'))

    
    # objective (eh.... how to create empty thing here)
    objective = edges[0].w * lp_x[0] 
    for i, edge in enumerate(edges):
        if i == 0:
            continue
        objective += edges[i].w * lp_x[i]
    prob += objective


    def delta_part(nodeset, edges):
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
        if len(nodeset) > 0 and len(nodeset) < len(pathtsp.graph.nodes()) and stcount in {0,2}:
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
def find_x_star(pathtsp, integer=False, b_cuts=None):
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



def find_B(pathtsp, x):
    """
    returns a list of cuts: all cuts in B(x)
    """
    B = set()
    # iterate over all subsets
    for nodeset in powerset(pathtsp.graph.nodes()):
        # only cuts with s in it and t not in it
        if pathtsp.s not in nodeset:
            continue
        if pathtsp.t in nodeset:
            continue

        # only cuts such that x(delta(nodeset)) < 3
        this_cut_value = 0
        edges = pathtsp.graph.edges()
        for i in delta(nodeset, edges):
            this_cut_value += x[edges[i]]
        
        if this_cut_value < 3:
            B.add(nodeset)
    return B



INFTY  = 100000000000

def get_y(pathtsp, B_family):
    # create a new family
    new_b = list(copy.deepcopy(B_family))
    new_b.append(tuple())
    final_b = tuple(pathtsp.graph.nodes())
    final_e = Edge(pathtsp.t, Node('dummydummy'), 0)
    value, y_star = get_y_recursive(pathtsp, new_b, final_b, final_e)
    return y_star

def get_y_recursive(pathtsp, B_family, B, e, dp_memo=None):
    """
    finds y recursively, by trying all possible values for where
    the previous 1-cut could be
    """
    print(f'e: {e}')
    if dp_memo is None:
        dp_memo = {}
    if (B, e) in dp_memo:
        return dp_memo[(B,e)]
    if len(B) == 0:
        # base case!!!
        return 0, {}
    edges = pathtsp.graph.edges() + [Edge(Node('dummy'),pathtsp.s,0)]
    minvalhere = INFTY
    minsolhere = None
    for B_prime in B_family:
        # check if B_prime is subset of B
        if not set(B_prime) < set(B):
            continue
        if len({e.n1,e.n2}.intersection(B_prime)) > 0:
            continue
        for ei in delta(B_prime, edges, emptyindex = len(edges)-1):
            e_prime = edges[ei]
            print(f'B: {B} B_prime: {B_prime} e: {e} e_prime: {e_prime}')
            # make sure we have a valid chain
            intersectgraph = set(B).difference(set(B_prime))
            if len({e_prime.n1,e_prime.n2}.intersection(intersectgraph)) == 0:
                continue
            # now call dp on this
            this_val, this_sol = get_y_recursive(pathtsp, B_family, B_prime, e_prime, dp_memo=dp_memo)
            # now solve the LP
            # we need to be careful here because we might have that the start node
            # is equal to the end node. that's fine tho.
            if {e_prime.n1,e_prime.n2}.intersection({e.n1,e.n2}) == intersectgraph:
                if len(intersectgraph) != 1:
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
                print(intersectgraph)
                print(e)
                print(B_prime)
                print(e_prime)
                (s,) = {e_prime.n1, e_prime.n2}.intersection(intersectgraph)
                (t,) = {e.n1, e.n2}.intersection(intersectgraph)
                sub_pathtsp = PathTSPProblem(induced_graph, s, t)
                print(induced_graph)
                print(pathtsp.graph)
                relevant_b_cuts = [set(b).intersection(intersectgraph) for b in B_family if set(b) < set(B) and set(B_prime) < set(b)]
                y_i = find_x_star(sub_pathtsp, b_cuts=relevant_b_cuts)

                val_y_i = sum([k.w*y_i[k] for k in y_i])

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

    dp_memo[(B, e)] = (minvalhere, copy.deepcopy(minsolhere))
    return dp_memo[(B, e)]




    



def integrality_gap_ptsp_graph(n):
    # look at page 2 of shmoys
    v = [Node('s'),Node('t')] + [Node(f't{i}') for i in range(n)] + [Node(f'b{i}') for i in range(n)]
    # t is top row b is bottom row
    e = [Edge(Node('s'), Node('t0'), 1), Edge(Node('s'),Node('b0'), 1), Edge(Node('t0'),Node('b0'), 1)] + \
        [Edge(Node('t'), Node(f't{n-1}'), 1), Edge(Node('t'),Node(f'b{n-1}'), 1), Edge(Node(f't{n-1}'),Node(f'b{n-1}'), 1)] + \
        [Edge(Node(f't{i}'), Node(f't{i+1}'), 1) for i in range(n-1)] + \
        [Edge(Node(f'b{i}'), Node(f'b{i+1}'), 1) for i in range(n-1)]
    s = SimpleGraph(v, e)
    
    s.to_dot('ssss.dot')

    m = s.to_metric()
    return m, e



