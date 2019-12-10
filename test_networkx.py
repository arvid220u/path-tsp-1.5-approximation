import networkx as nx
import matplotlib.pyplot as plt
import pydot

G=nx.Graph()

G.add_node(1)

G.add_nodes_from([2,3])

G.add_edge(1,2)

print(G.edges())

# nx.draw(G)

ba=nx.barabasi_albert_graph(100,5)
# nx.draw_circular(ba)

nx.nx_pydot.write_dot(G,'file.dot')

plt.show()
