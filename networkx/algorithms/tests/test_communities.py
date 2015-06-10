#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestCommunity:

    def test_louvain_undirected(self):
        G = nx.Graph()
        G.add_edges_from([(0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 4),
                          (1, 7), (2, 4), (2, 5), (2, 6), (3, 7), (4, 10),
                          (5, 7), (5, 11), (6, 7), (6, 11), (8, 9), (8, 10),
                          (8, 11), (8, 14), (8, 15), (9, 12), (9, 14),
                          (10, 11), (10, 12), (10, 13), (10, 14), (11, 13)])
        p = nx.louvain(G)

    @raises(nx.NetworkXError)
    def test_louvain_directed(self):
        G = nx.DiGraph()
        G.add_weighted_edges_from([('0', '3', 3), ('0', '1', -5),
                                   ('1', '0', -2), ('0', '2', 2),
                                   ('1', '2', -3), ('2', '3', 1)])
        assert_raises(nx.NetworkXError, nx.louvain(G))

    @raises(nx.NetworkXError)
    def test_louvain_no_edges(self):
        G = nx.Graph()
        G.add_nodes_from([1, 2, 3, 4])
        assert_raises(nx.NetworkXError, nx.louvain(G))
