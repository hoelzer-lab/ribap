#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ILP file of RIBAP. This file contains several classes that are needed
for the generation of pairwise ILPs.

Usage:
ILP.py [options] <tsv_file>

Options:
    -h, --help                                      Show this neat help message and exit.
    -s SELFWEIGHT, --self SELFWIGHT                 Weight of a self-edge in the matching. [default: 0.0]
    -a ALPHA, --alpha ALPHA                         Alpha coefficient for the rearrangement/sequence balance. [default: 0.5]
    -x TOLERANCE, --fix_adj_tolerance TOLERANCE     Weight tolerance for fixing gene adjacencies. [default: 0.0]
    -max, --maximal                                 Force maximal matching.
    -ms, --matching_size                            Matching size in the F.O.
    -i, --indel                                     Use measurement of InDel events.

"""

__author__ = "Kevin Lamkiewicz"

from itertools import chain
import sys
import re

import networkx as nx
from networkx.algorithms import bipartite
from networkx import connected_components

from docopt import docopt



class GSGraph(nx.Graph):

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
            return "A%d" % idx
        else:
            return "B%d" % idx


class FFProblem():

    def __init__(self, filename):
        self.G = None
        self.n1 = 0
        self.n2 = 0
        self._read_ff_file(filename, isolated=True)

    def size(self, genome):
        return self.n1 if genome == 0 else self.n2

    def g1_degrees_sorted(self):
        g1_nodes = set(n for n, d in self.G.nodes(data=True) if d['bipartite'] == 0)
        return sorted(self.G.degree_iter(g1_nodes), key=itemgetter(1))

    def _read_ff_file(self, filename, isolated=True):
        G = GSGraph()
        last = 0
        max_g2 = 0
        g1 = 0  # will store the max g1 at end
        for l in open(filename):
            if not l.strip():
                continue
            g1, g2, orient, w = l.split()
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

class ILPGenerator():

    def __init__(self, tsvFile):
        
        self.dirname = '/'.join(tsvFile.split('/')[:-2])
        self.dirname = self.dirname + "/ilp/" if self.dirname else "ilp/"
        self.basename = tsvFile.split('/')[-1].split('.')[0]
        self.out = self.dirname+self.basename 
        
        
        
    def atoi(self, text):
        return int(text) if text.isdigit() else text


    def natural_keys(self, text):
        return [self.atoi(c) for c in re.split('(\d+)', text)]


    def all_extremities(self, P):
        for g in range(2):
            for i in range(1, P.size(g) + 1):
                yield GSGraph.node_name(i, g) + "h"
                yield GSGraph.node_name(i, g) + "t"


    def all_extremities_of_genome(self, P, g):
        for i in range(1, P.size(g) + 1):
            yield GSGraph.node_name(i, g) + "h"
            yield GSGraph.node_name(i, g) + "t"


    def weighted_edge_var(self, u, v, w):
        if w > 0:
            return [self.edge_var(u, 't', v, 't'), self.edge_var(u, 'h', v, 'h')]
        else:
            return [self.edge_var(u, 't', v, 'h'), self.edge_var(u, 'h', v, 't')]


    def edge_var(self, u, ext_u, v, ext_v):
        #return "x_%s_%s" % tuple(sorted([u + ext_u, v + ext_v]))
        sortedTuple = tuple(sorted([u + ext_u, v + ext_v]))
        return f"x_{sortedTuple[0]}_{sortedTuple[1]}"


    def test_edge(self, G, a, b, positive, nodes, tolerance=1):
        v_a = GSGraph.node_name(a, 0)
        v_b = GSGraph.node_name(b, 1)
        if v_a in nodes or v_b in nodes or not G.has_edge(v_a, v_b):
            return False
            #print >> sys.stderr,v_a,v_b, "Edge"
        w = G[v_a][v_b]['weight']
        #print >> sys.stderr, w, positive
        if (w > 0 and positive) or (w < 0 and not positive):  # either w>0 and positive == True, or w<0 and positive == False
            # test all neightbours
            for u in G.neighbors(v_a):
                w_u = G[v_a][u]['weight']
                if (w_u > 0 and positive or w_u < 0 and not positive) and abs(w_u) > abs(w) * tolerance:
                    return False
            for u in G.neighbors(v_b):
                w_u = G[v_b][u]['weight']
                if (w_u > 0 and positive or w_u < 0 and not positive) and abs(w_u) > abs(w) * tolerance:
                    return False
        return True


    def fix_adjacencies(self, G, edges, nodes, tolerance=1):
        new_edges = []
        # print >> sys.stderr, "EDGES:", edges
        # print >> sys.stderr, "NODES:", nodes
        while len(edges) > 0:
            # pick one edge:
            a, b = edges.pop()
            # edge vertices are not always in A,B order, fix:
            if a[0] == 'B':
                b, a = a, b
                # get weight:
            w = G[a][b]['weight']
            a, b = map(int, [a[1:], b[1:]])
            if w > 0:
                adj = [(a + 1, b + 1, True), (a - 1, b - 1, True)]
            else:
                adj = [(a - 1, b + 1, False), (a + 1, b - 1, False)]
            for va, vb, positive in adj:
                r = test_edge(G, va, vb, positive, nodes, tolerance)
                if r:
                    v_a = GSGraph.node_name(va, 0)
                    v_b = GSGraph.node_name(vb, 1)
                    # print >> sys.stderr, "add", v_a, v_b
                    # remove edges:
                    for u in G.neighbors(v_a):
                        if u != v_b:
                            G.remove_edge(v_a, u)
                    for u in G.neighbors(v_b):
                        if u != v_a:
                            G.remove_edge(v_b, u)

                    new_edges.append([v_a, v_b])
                    edges.append([v_a, v_b])
                    nodes.append(v_a)
                    nodes.append(v_b)

        #print >> sys.stderr, "Fixed %d edges." % len(new_edges)
        return new_edges

    def get_stream(self, u, components, streams):
        for idx, component in components.items():
            if u in component:
                return streams[idx]

    def broadcast_streams(self, streams, message):
        for idx, stream in streams.items():
            stream.write(message)

    def generate_lp(self, P, alpha=0.5, self_edge_cost=0, MATCHING_SIZE=False, MAXIMAL_MATCHING=False, tolerance=0, INDEL=False):
        G = P.G
        outputStreams = {}
        disjointComponents = {idx:c for idx,c in enumerate(connected_components(G)) if len(c)}
        #component2Stream = {c : idx for idx,c in disjointComponents.items()}
        #exit(0)
        for idx, _ in disjointComponents.items():
            stream = open(f"{self.out}_{idx}.ilp",'w')
            outputStreams[idx] = stream

        self.broadcast_streams(outputStreams, "Minimize\nobj: ")    
        if not MATCHING_SIZE:
            self.broadcast_streams(outputStreams, f"{(P.n1 + P.n2) / 2.0} g ")
            
        
        # Find K2 components, to fix edges:
        k2 = [c for c in connected_components(G) if len(c) == 2]

        # try to fix adjacencies:
        if tolerance > 0:
            #print >> sys.stderr, "Fixing edges..."
            new_edges = fix_adjacencies(G, list(k2), list(chain.from_iterable(k2)), tolerance)
            k2 += new_edges

        k2_vertices = list(chain.from_iterable(k2))
        constraintID = 1
        # objective function:
        #print >> sys.stderr, "Printing ILP ..."
        #print("Minimize")
        #print("obj:", end=' ')
        #if not MATCHING_SIZE:
            # print P.n1 + P.n2, "g",
            #print(P.n1 + P.n2) / 2.0, "g",  # new alpha
        #    print(f"{(P.n1 + P.n2) / 2.0} g", end=" ")  # new alpha
            #print(f"{(P.n1 + P.n2) / 2.0}", end=" ")  # new alpha

            # Since the weight of "parallel edges" is the same, and I need to divide by 2, and also parallel edges are
        # forced to be together, I only include one of each in the OF (only edges[0] below)
        for u, v in G.edges():
            stream = self.get_stream(u,disjointComponents,outputStreams)
            w = G[u][v]['weight']
            edges = self.weighted_edge_var(u, v, w)
            # edge_w = 2 - abs(w) if MATCHING_SIZE else -abs(w)
            edge_w = (1 if MATCHING_SIZE else 0) + (alpha - 1) * abs(w)
            stream.write(f"{edge_w} {edges[0]} ")
            #print("%+f %s" % (edge_w, edges[0]), end=" ")
            #print(f"{edge_w} {edges[0]}", end=" ")
        
        # 28.06.2016, Kevin:
        # Changed the self_edge_cost to alpha in order to implement the modification proposed
        # in Braga et al. 2011 - DCJ with Indels
        if INDEL:
            self_edge_cost = alpha

        if self_edge_cost != 0:
            for g in range(2):  # each genome
                for i in range(1, P.size(g) + 1):
                    u = GSGraph.node_name(i, g)
                    # ignore k2 vertices, they are already fixed;
                    if u in k2_vertices:
                        continue
                    self.get_stream(u,disjointComponents,outputStreams)
                    #print(f"{self_edge_cost} {edge_var(u, "t", u, "h")}", end=" ")
                    stream.write(f"+ {self_edge_cost} {self.edge_var(u, 't', u, 'h')} ")
                    
                    #print("%+.2f %s" % (self_edge_cost, self.edge_var(u, "t", u, "h")), end =" ")                    

        # print "-", " - ".join(["%.2f z_%s" % (cycle_weight, x) for x in
        # all_extremities_of_genome(P, 0) if x[-1] == 'h']),
        #print("-", " - ".join(["%.3f z_%s" % (alpha, x) for x in self.all_extremities_of_genome(P, 0) if x[-1] == 'h']), end=" ")
        for x in self.all_extremities_of_genome(P,0):
            if x[-1] == 'h':
                stream = self.get_stream(x[:-1], disjointComponents, outputStreams)
                stream.write(f"- {alpha} z_{x} ")
        
        if INDEL:
            for x in sorted(G.nodes()):
                if x not in k2_vertices:
                    stream = self.get_stream(x, disjointComponents, outputStreams)
                    stream.write(f"+ {alpha} b_{x} ")
            #print("+", " + ".join("%.3f b_%s" % (alpha, x) for x in sorted(G.nodes()) if x not in k2_vertices), end=" ")
        #print()

        # constraints:
        self.broadcast_streams(outputStreams, "\nSubject To\n")
        #for idx, stream in outputStreams.items():
            #print("Subject To")
        #    stream.write("\nSubject To\n")

        if not MATCHING_SIZE:
            #print(f"c{constraintID}: g = 1")
            self.broadcast_streams(outputStreams, f"c{constraintID}: g = 1\n")
            #stream.write(f"c{constraintID} g = 1\n")
            constraintID += 1

            # Consistency:
            #print("\ consistency (parallel edges)")
            stream.write("\ consistency (parallel edges)\n")
        

        for v, u in G.edges():
            w = G[v][u]['weight']
            stream = self.get_stream(u, disjointComponents, outputStreams)
            if u in k2_vertices and v in k2_vertices:
                for e in self.weighted_edge_var(u, v, w):
                    #print e, "= 1")
                    #print(f"c{constraintID}: {e} = 1")
                    stream.write(f"c{constraintID}: {e} = 1\n")
                    constraintID += 1
            else:
                #print("c"+str(constraintID) + ": " + " - ".join(self.weighted_edge_var(u, v, w)), "= 0")
                stream.write("c"+str(constraintID) + ": " + " - ".join(self.weighted_edge_var(u, v, w)) + "= 0\n")
                constraintID += 1
        
        # One edge per vertex:
        self.broadcast_streams(outputStreams, "\ degree 1 for each vertex\n")
        #for idx, stream in outputStreams.items():
        #    stream.write("\ degree 1 for each vertex\n")
        #print ("\ degree 1 for each vertex")
        
        for g in range(2):  # each genome
            for i in range(1, P.size(g) + 1):
                u = GSGraph.node_name(i, g)
                # ignore k2 vertices, they are already fixed;
                if u in k2_vertices:
                    continue
                stream = self.get_stream(u, disjointComponents, outputStreams)
                nodes_h = []
                nodes_t = []
                for v in G.neighbors(u):
                    edges = self.weighted_edge_var(u, v, G[v][u]['weight'])
                    nodes_t.append(edges[0])
                    nodes_h.append(edges[1])
                    # self edges: (edges to complete the matching for unsaturated vertices
                nodes_t.append(self.edge_var(u, "t", u, "h"))
                nodes_h.append(self.edge_var(u, "t", u, "h"))
                # print "# %st" % u
                #print(" + ".join(nodes_t), " = 1", "    \\", u + "t")
                #print("c"+str(constraintID) + ": "+" + ".join(nodes_t), " = 1")
                stream.write("c"+str(constraintID) + ": "+" + ".join(nodes_t) + " = 1\n")
                constraintID += 1
                # print "# %sh" % u
                #print("c"+str(constraintID) + ": "+" + ".join(nodes_h), " = 1")
                stream.write("c"+str(constraintID) + ": "+" + ".join(nodes_h) + " = 1\n")
                constraintID += 1

        # Maximal matching:
        if MAXIMAL_MATCHING:
            self.broadcast_streams(outputStreams, "\ Maximal matching\n")
            #print("\ Maximal matching")
            for u, v in G.edges():
                if u in k2_vertices:
                    continue
                stream = self.get_stream(u, disjointComponents, outputStreams)
                #print("c"+str(constraintID) + ": " + "%s + %s <= 1" % (self.edge_var(u, "t", u, "h"), self.edge_var(v, "t", v, "h")))
                stream.write("c"+str(constraintID) + ": " + "%s + %s <= 1\n" % (self.edge_var(u, "t", u, "h"), self.edge_var(v, "t", v, "h")))
                constraintID += 1
        
        # Adjacent vertices, same label:
        m = (P.n1 + P.n2) * 2
        #print("\ Label - adj vertices with same label")
        self.broadcast_streams(outputStreams, "\ Label - adj vertices with same label\n")
        # matching edges:
        for u, v in G.edges():
            stream = self.get_stream(u, disjointComponents, outputStreams)
            edges = self.weighted_edge_var(u, v, G[u][v]['weight'])
            for i in [0, 1]:
                x, v1, v2 = edges[i].split("_")
                # print
                if u in k2_vertices:
                    #print(f"c{constraintID}: y_{v1} - y_{v2} = 0")
                    stream.write(f"c{constraintID}: y_{v1} - y_{v2} = 0\n")
                    constraintID += 1
                else:
                    #print("y_%s - y_%s  + %d %s <= %d" % (v1, v2, m, edges[i], m))
                    #print(f"c{constraintID}: y_{v1} - y_{v2}  + {m} {edges[i]} <= {m}")
                    stream.write(f"c{constraintID}: y_{v1} - y_{v2}  + {m} {edges[i]} <= {m}\n")
                    constraintID += 1
                    #print("y_%s - y_%s  + %d %s <= %d" % (v2, v1, m, edges[i], m))
                    #print(f"c{constraintID}: y_{v2} - y_{v1}  + {m} {edges[i]} <= {m}")
                    stream.write(f"c{constraintID}: y_{v2} - y_{v1}  + {m} {edges[i]} <= {m}\n")
                    constraintID += 1
                    # self edges:
        
        for g in range(2):
            for i in range(1, P.size(g) + 1):
                u = GSGraph.node_name(i, g)
                # k2 vertices do not have self edges:
                if u in k2_vertices:
                    continue
                stream = self.get_stream(u, disjointComponents, outputStreams)
                #print("y_%sh - y_%st  + %d %s <= %d" % (u, u, m, self.edge_var(u, 'h', u, 't'), m), "\  Self edge for ", u)
                #print("y_%sh - y_%st  + %d %s <= %d" % (u, u, m, self.edge_var(u, 'h', u, 't'), m), "\  Self edge for ", u)
                #print(f"c{constraintID}: y_{u}h - y_{u}t  + {m} {self.edge_var(u, 'h', u, 't')} <= {m}")
                stream.write(f"c{constraintID}: y_{u}h - y_{u}t  + {m} {self.edge_var(u, 'h', u, 't')} <= {m}\n")
                constraintID += 1
                #print("y_%st - y_%sh  + %d %s <= %d" % (u, u, m, self.edge_var(u, 'h', u, 't'), m))
                #print(f"c{constraintID}: y_{u}t - y_{u}h  + {m} {self.edge_var(u, 'h', u, 't')} <= {m}")
                stream.write(f"c{constraintID}: y_{u}t - y_{u}h  + {m} {self.edge_var(u, 'h', u, 't')} <= {m}\n")
                constraintID += 1
                #        nodes_t.append(edge_var(u, "t", u, "h"))

        # adjacency edges:
        #print()
        for g in range(2):
            for i in range(1, P.size(g)):
                #print("y_%s - y_%s = 0" % (GSGraph.node_name(i, g) + "h", GSGraph.node_name(i + 1, g) + "t"))
                #print("c"+str(constraintID)+": " + "y_%s - y_%s = 0" % (GSGraph.node_name(i, g) + "h", GSGraph.node_name(i + 1, g) + "t"))
                self.broadcast_streams(outputStreams, "c"+str(constraintID)+": " + "y_%s - y_%s = 0\n" % (GSGraph.node_name(i, g) + "h", GSGraph.node_name(i + 1, g) + "t"))
                constraintID += 1
                # Assuming circular genome:
            #print("y_%s - y_%s = 0" % (GSGraph.node_name(P.size(g), g) + "h", GSGraph.node_name(1, g) + "t"))
            self.broadcast_streams(outputStreams, "c"+str(constraintID)+": " + "y_%s - y_%s = 0\n" % (GSGraph.node_name(P.size(g), g) + "h", GSGraph.node_name(1, g) + "t"))
            #print("c"+str(constraintID)+": " + "y_%s - y_%s = 0" % (GSGraph.node_name(P.size(g), g) + "h", GSGraph.node_name(1, g) + "t"))
            constraintID += 1

        # cycle counter:
        for c, e in enumerate(self.all_extremities_of_genome(P, 0)):
            if e[-1] == 't':
                continue
            self.broadcast_streams(outputStreams, "c"+str(constraintID)+": " + "%d z_%s - y_%s <= 0\n" % (c + 1, e, e))
            #print("c"+str(constraintID)+": " + "%d z_%s - y_%s <= 0" % (c + 1, e, e))
            constraintID += 1

        ##################
        # Second Approach: We try to measure the contigousness of selfedges
        ##################
        if INDEL:
            for g in range(2):
                for i in range(1, P.size(g) + 1):
                    currentGene = GSGraph.node_name(i, g)
                    if not (currentGene in k2_vertices):
                        if i == P.size(g):
                            # continue
                            nextGene = GSGraph.node_name(1, g)
                        else:
                            nextGene = GSGraph.node_name(i + 1, g)
                        if not nextGene in k2_vertices:
                            self.broadcast_streams(outputStreams, "c"+str(constraintID)+": " + "%s - %s - b_%s <= 0\n" % (self.edge_var(currentGene, 'h', currentGene, 't'), self.edge_var(nextGene, 'h', nextGene, 't'), currentGene))
                            #print("c"+str(constraintID)+": " + "%s - %s - b_%s <= 0" % (self.edge_var(currentGene, 'h', currentGene, 't'), self.edge_var(nextGene, 'h', nextGene, 't'), currentGene))
                            constraintID += 1
                        else:
                            self.broadcast_streams(outputStreams, "c"+str(constraintID)+": " + "%s - b_%s = 0\n" % (self.edge_var(currentGene, 'h', currentGene, 't'), currentGene))
                            #print("c"+str(constraintID)+": " + "%s - b_%s = 0" % (self.edge_var(currentGene, 'h', currentGene, 't'), currentGene))
                            constraintID += 1
                

        # Bounds:
        self.broadcast_streams(outputStreams, "\nBounds\n")
        #print()
        #print("Bounds")
        # Node labels:
        c = 1
        for e in self.all_extremities(P):
            self.broadcast_streams(outputStreams, "y_%s <= %d\n" % (e, c))
            #print("y_%s <= %d" % (e, c))
            c += 1
        #print()
        # Variable types:
        #print("Binary")
        self.broadcast_streams(outputStreams, "\nBinary\n")
        # matching edges
        for u, v in G.edges():
            stream = self.get_stream(u, disjointComponents, outputStreams)
            w = G[u][v]['weight']
            edges = self.weighted_edge_var(u, v, w)
            #print("\n".join(edges))
            stream.write("\n".join(edges)+"\n")
            # self edges:

        for g in range(2):
            for i in range(1, P.size(g) + 1):
                u = GSGraph.node_name(i, g)
                # k2 vertices do not have self edges:
                if u in k2_vertices:
                    continue
                stream = self.get_stream(u, disjointComponents, outputStreams)
                #print(self.edge_var(u, 'h', u, 't'))
                stream.write(self.edge_var(u, 'h', u, 't')+"\n")
            # z_ : cycle counters (only vertices in A of type 'head')

        self.broadcast_streams(outputStreams, "\n".join("z_" + x for x in self.all_extremities_of_genome(P, 0) if x[-1] == 'h'))
        #print("\n".join("z_" + x for x in self.all_extremities_of_genome(P, 0) if x[-1] == 'h'))
        #print()
        if INDEL:
            self.broadcast_streams(outputStreams, "\n"+"\n".join("b_" + x for x in sorted(G.nodes()) if x not in k2_vertices))
            #print("\n".join("b_" + x for x in sorted(G.nodes()) if x not in k2_vertices))
        #if not MATCHING_SIZE:
        #    print("g")
        self.broadcast_streams(outputStreams, "\nGeneral\n"+"\n".join("y_" + x for x in self.all_extremities(P))+"\nEnd\n")
        #print("General")
        #print("\n".join("y_" + x for x in self.all_extremities(P)))
        #print
        # print "\n".join("k_" + x for x in sorted(GSGraph.nodes(G)) if x not in k2_vertices)
        #print("End")


def parse_arguments(param):
    s = float(param['--self'])
    alpha = float(param['--alpha'])
    matching_size = param['--matching_size']
    maximal = param['--maximal']
    fix_adj_tolerance = float(param['--fix_adj_tolerance'])
    indel = param['--indel']

    tsvFile=param['<tsv_file>']

    return [tsvFile, s, alpha, matching_size, maximal, fix_adj_tolerance, indel]

if __name__ == '__main__':


    param = docopt(__doc__)
    tsvFile, s, alpha, matching_size, maximal, fix_adj_tolerance, indel = parse_arguments(param)
    #print >> sys.stderr, param

    # Read problem
    P = FFProblem(tsvFile)

    # Filter edges:
    #print >> sys.stderr, "Filtering edges..."
    #if param.f > 0:
    #    ff_graph.filter_edges_from_max_edge(P.G, param.f)

    # print ILP:
    #print >> sys.stderr, "Generating ILP..."
    ilpGen = ILPGenerator(tsvFile)
    ilpGen.generate_lp(P, self_edge_cost=s, alpha=alpha, MATCHING_SIZE=matching_size,
                MAXIMAL_MATCHING=maximal, tolerance=fix_adj_tolerance, INDEL=indel)