import networkx as nx

class Node:
    def __init__(self, label):
        self.label = label
    
    def __hash__(self):
        return hash(self.label)

    def __eq__(self, other):
        return self.label == other.label

    def __ne__(self, other):
        return not(self == other)
    
    def __str__(self):
        return str(self.label)

    def __lt__(self,other):
        return self.label < other.label
    
    def __le__(self,other):
        return self.label <= other.label
    
    def __repr__(self):
        return self.__str__()

class Edge:
    def __init__(self, n1, n2, w):
        self.n1 = n1
        self.n2 = n2
        self.w = w
    
    def __str__(self):
        return str(self.n1) + "->" + str(self.n2) + "(" + str(self.w) + ")"

    def __hash__(self):
        return hash((self.n1, self.n2, self.w))

    def __eq__(self, other):
        return (self.n1, self.n2, self.w) == (other.n1, other.n2, other.w)

    def __ne__(self, other):
        return not(self == other)
    
    def __repr__(self):
        return self.__str__()

class SimpleGraph:
    def __init__(self, nodes, edge_list):
        self.nodes = nodes
        self.edge_list = []
        self.infty = 10000000
        for e in edge_list:
            self.edge_list.append(Edge(e.n1,e.n2,e.w))
            self.edge_list.append(Edge(e.n2,e.n1,e.w))
    
    def to_metric(self):
        m = MetricGraph(self.nodes)
        # solve all pairs shortest path problem
        # floyd warshall
        dist = {nn: {n: self.infty for n in self.nodes} for nn in self.nodes}
        for e in self.edge_list:
            dist[e.n1][e.n2] = e.w
        for n in self.nodes:
            dist[n][n] = 0
        for k in self.nodes:
            for i in self.nodes:
                for j in self.nodes:
                    if dist[i][j] > dist[i][k] + dist[k][j]:
                        dist[i][j] = dist[i][k] + dist[k][j]
        # now create edges for every distance
        for n1 in self.nodes:
            for n2 in self.nodes:
                m.set_distance(n1,n2,dist[n1][n2])
        m.assert_metric()
        return m
    
    def to_dot(self, path):
        # create networkx
        G = nx.Graph()
        for n in self.nodes:
            G.add_node(n)
        for e in self.edge_list:
            G.add_edge(e.n1, e.n2)
        nx.nx_pydot.write_dot(G, path)


class MetricGraph:
    """
    this graph can for example be created from a shortest path
    run in another graph
    """
    def __init__(self, nodes):
        """
        assumes nodes are from 0 to n
        """
        self.adj_matrix = {nn: {n: None for n in nodes} for nn in nodes}
        for n in nodes:
            self.adj_matrix[n][n] = Edge(n,n,0)
    
    def distance(self, x, y):
        return self.adj_matrix[x][y].w
    
    def set_distance(self, x, y, d):
        self.adj_matrix[x][y] = Edge(x,y,d)
    
    def nodes(self):
        return sorted(list(self.adj_matrix.keys()))
    
    def edges(self):
        l = []
        for n1 in self.nodes():
            for n2 in self.nodes():
                if n1 < n2:
                    l.append(self.adj_matrix[n1][n2])
        return l
    
    def assert_metric(self):
        """
        assserts that the graph is indeed metric
        """
        for x in self.adj_matrix:
            for y in self.adj_matrix:
                if type(self.adj_matrix[x][y]) != type(Edge(0,0,0)):
                    raise RuntimeError("not all edges are set!!")

        for x in self.adj_matrix:
            for y in self.adj_matrix:
                for z in self.adj_matrix:
                    if self.distance(x,y) + self.distance(y,z) < self.distance(x,z):
                        raise RuntimeError("violate triangle ineq!!")
        
        for x in self.adj_matrix:
            for y in self.adj_matrix:
                if self.distance(x,y) != self.distance(x,y):
                    raise RuntimeError("violates symmetric property!!")

        for x in self.adj_matrix:
            for y in self.adj_matrix:
                if x == y:
                    if self.distance(x,y) != 0:
                        raise RuntimeError("distance to itself must be 0!!")
                else:
                    if self.distance(x,y) <= 0:
                        raise RuntimeError("for x != y distance needs to be positive")

    def __str__(self):
        # print self as adjacency matrix
        s = ""
        for n1 in self.nodes():
            for n2 in self.nodes():
                s += str(self.adj_matrix[n1][n2]) + " "
            s += '\n'
        return s
    

    def subset_graph(self, v):
        """
        returns the induced graph on the subset of vertices v
        """
        induced_graph = MetricGraph(v)
        for v1 in v:
            for v2 in v:
                induced_graph.set_distance(v1, v2, self.distance(v1, v2))
        return induced_graph
