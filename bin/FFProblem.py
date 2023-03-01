from GSGraph import GSGraph

class FFProblem():
  
    def __init__(self, strains, sims):
        self.G = None
        self.n1 = 0
        self.n2 = 0
        self._read_ff_file(strains, sims, isolated=True)

    def size(self, genome):
        return self.n1 if genome == 0 else self.n2

    def g1_degrees_sorted(self):
        g1_nodes = set(n for n, d in self.G.nodes(data=True) if d['bipartite'] == 0)
        return sorted(self.G.degree_iter(g1_nodes), key=itemgetter(1))

    def _read_ff_file(self, strains, sims, isolated=True):
        G = GSGraph()
        last = 0
        max_g2 = 0
        g1 = 0  # will store the max g1 at end
        for gene_sim in sims:
        #for l in open(filename):
            #if not l.strip():
                #continue
            g1, g2, orient, w = gene_sim
            g1 = int(g1)
            g2 = int(g2)
            w = float(w) * int(orient)
            if isolated:
                if g1 > last:
                    for i in range(last + 1, g1 + 1):
                        G.add_genome_node(i, 0)
                    last = g1
            G.add_genome_edge(g1, g2, w)
            if g2 > max_g2:
                max_g2 = g2
        if isolated:
            for i in range(1, max_g2 + 1):
                G.add_genome_node(i, 1)
        self.G = G
        self.n1 = g1
        self.n2 = max_g2