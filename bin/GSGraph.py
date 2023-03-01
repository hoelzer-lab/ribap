import networkx as nx

class GSGraph(nx.Graph):
    """
    
    """

    def add_genome_node(self, idx, genome):
        self.add_node(GSGraph.node_name(idx, genome), bipartite=genome)

    def add_genome_edge(self, g1, g2, w):
        self.add_genome_node(g1, 0)
        self.add_genome_node(g2, 1)
        self.add_edge(GSGraph.node_name(g1, 0), GSGraph.node_name(g2, 1), weight=w)

    @staticmethod
    def next_node(u):
        return u[0] + str(int(u[1:]) + 1)

    @staticmethod
    def previous_node(u):
        return u[0] + str(int(u[1:]) - 1)

    @staticmethod
    def node_name(idx, genome):
        if genome == 0:
            return f"A{idx}"
        else:
            return f"B{idx}"