from path_tsp_approx import *


def try_lp():
    # this works! p cool.
    x1 = pulp.LpVariable('x1',lowBound=0)
    x2 = pulp.LpVariable('x2',lowBound=0)
    x3 = pulp.LpVariable('x3',lowBound=0)
    prob = pulp.LpProblem('testproblem',sense=pulp.LpMinimize)
    prob += x1 + x2 >= 1
    constraint = x1 + 2*x2
    constraint = constraint == 3
    prob += constraint

    objective = 2*x3
    objective += x2 + x3
    prob += objective
    status = prob.solve()

    print(pulp.LpStatus[status])
    print(pulp.value(x1))
    print(pulp.value(objective))
    print(prob)

def try_hk():
    # create a graph.
    v = [Node(0),Node(1),Node(2),Node(3)]
    e = [Edge(v[0],v[1],1),
         Edge(v[1],v[2],1),
         Edge(v[2],v[3],1),
         Edge(v[0],v[3],1)]
    s = SimpleGraph(v, e)
    m = s.to_metric()
    print(m)
    p = PathTSPProblem(m, v[0],v[3])
    hk = held_karp_lp(p)
    print(hk)

def try_x_star():
    # create a graph.
    v = [Node(0),Node(1),Node(2),Node(3)]
    e = [Edge(v[0],v[1],1),
         Edge(v[1],v[2],1),
         Edge(v[2],v[3],1),
         Edge(v[0],v[3],1)]
    s = SimpleGraph(v, e)
    m = s.to_metric()
    print(m)
    p = PathTSPProblem(m, v[0],v[3])
    find_x_star(p)

def try_x_star_hard():
    # create a graph.
    m, e = integrality_gap_ptsp_graph(4)
    p = PathTSPProblem(m, Node('s'), Node('t'))
    x_star = find_x_star(p, integer=True)
    x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    print(x_star_supp)
    # s = SimpleGraph(m.nodes(), list(x_star_supp.keys()) + e)
    # s.to_dot('ssss.dot')

def try_b_x_star():
    m, e = integrality_gap_ptsp_graph(4)
    p = PathTSPProblem(m, Node('s'), Node('t'))
    x_star = find_x_star(p)
    # print(x_star)
    B = find_B(p, x_star)
    x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    print(x_star_supp)
    print(B)

def try_get_y():
    m, e = integrality_gap_ptsp_graph(3)
    p = PathTSPProblem(m, Node('s'), Node('t'))
    x_star = find_x_star(p)
    # print(x_star)
    B = find_B(p, x_star)
    x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    print(x_star_supp)
    print(B)
    y_star = get_y(p, B)
    print(y_star)
    # this is so interesting!!! it find the optimum solution :0000

def try_get_y_2():
    m, e = y_not_0_graph(3)
    p = PathTSPProblem(m, Node('s'), Node('t'))
    x_star = find_x_star(p)
    # print(x_star)
    B = find_B(p, x_star)
    x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    print(x_star_supp)
    print(B)
    y_star = get_y(p, B)
    print(y_star)

def try_get_y_3():
    m, e = integrality_gap_ptsp_graph_perturbed(3)
    p = PathTSPProblem(m, Node('s'), Node('t'))
    x_star = find_x_star(p)
    # print(x_star)
    B = find_B(p, x_star)
    x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    print(x_star_supp)
    print(B)
    y_star = get_y(p, B)
    print(y_star)
    # this is so interesting!!! it find the optimum solution :0000

def try_k_branch_graph():
    m, e = k_branch_graph(4, 4)
    p = PathTSPProblem(m, Node('s'), Node('t'))
    x_star = find_x_star(p)
    print(x_star)
    # B = find_B(p, x_star)
    # x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    # print(x_star_supp)
    # print(B)
    # y_star = get_y(p, B)
    # print(y_star)
    # this is so interesting!!! it find the optimum solution :0000

def find_B_of_k_branch():
    import json
    with open('t.json','r') as f:
        x_star_str = json.load(f)
    # convert x_star_str into actual edges
    x_star = {}
    dd = None
    for k in x_star_str:
        n1, p2 = k.split('->')
        n2, w = p2.rstrip(')').split('(')
        w = int(w)
        e = Edge(Node(n1),Node(n2),w)
        x_star[e] = x_star_str[k]
        dd = e
    
    m, e = k_branch_graph(4, 4)
    p = PathTSPProblem(m, Node('s'), Node('t'))

    B = find_B(p, x_star)
    x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    print(x_star_supp)
    print(B)

def try_final():
    n = 3
    m, e = final_try_graph(n)
    p = PathTSPProblem(m, Node('0'), Node(f'{3*n}'))
    x_star = find_x_star(p)
    print(x_star)
    # B = find_B(p, x_star)
    x_star_supp = {d: x_star[d] for d in x_star if x_star[d] != 0}
    print(x_star_supp)
    # print(B)
    # y_star = get_y(p, B)
    # print(y_star)


if __name__ == '__main__':
    # try_lp()
    # try_hk()
    # try_x_star_hard()
    # try_b_x_star()
    # try_get_y_3()
    # try_get_y_2()
    # a = Node(0)
    # b = Node(1)
    # g = MetricGraph([a,b])
    # print(g.nodes())
    # print(a in g.nodes())
    # try_k_branch_graph()
    # find_B_of_k_branch()
    try_final()
    pass
    