import networkx as nx
import matplotlib.pyplot as plt


G = nx.DiGraph()
G.add_node("A")
G.add_node("B")
G.add_node("C")
G.add_node("D")
G.add_node("E")
G.add_node("F")
G.add_node("G")
G.add_edge("A","B")
G.add_edge("B","C")
G.add_edge("C","E")
G.add_edge("C","F")
G.add_edge("D","E")
G.add_edge("F","G")

print(G.nodes())
print(G.edges())

pos = nx.spring_layout(G)

nx.draw_networkx_nodes(G, pos)
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, edge_color='r', arrows = True)

plt.show()