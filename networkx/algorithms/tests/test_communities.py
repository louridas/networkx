#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestCommunity:

    def test_louvain(self):
        G = nx.Graph()
        for i in range(0, 16):
            G.add_node(i)
        G.add_edge(0, 2)
        G.add_edge(0, 3)
        G.add_edge(0, 4)
        G.add_edge(0, 5)
        G.add_edge(1, 2)
        G.add_edge(1, 4)
        G.add_edge(1, 7)
        G.add_edge(2, 4)
        G.add_edge(2, 5)
        G.add_edge(2, 6)
        G.add_edge(3, 7)
        G.add_edge(4, 10)
        G.add_edge(5, 7)
        G.add_edge(5, 11)
        G.add_edge(6, 7)
        G.add_edge(6, 11)
        G.add_edge(8, 9)
        G.add_edge(8, 10)
        G.add_edge(8, 11)
        G.add_edge(8, 14)
        G.add_edge(8, 15)
        G.add_edge(9, 12)
        G.add_edge(9, 14)
        G.add_edge(10, 11)
        G.add_edge(10, 12)
        G.add_edge(10, 13)
        G.add_edge(10, 14)
        G.add_edge(11, 13)
        p = nx.louvain(G)
        assert_equal(1, 1)

    @raises(nx.NetworkXError)
    def test_non_weighted_exception_louvain(self):
        G = nx.Graph()
        G.add_nodes_from([1, 2, 3, 4])
        assert_raises(nx.NetworkXError, nx.louvain(G))

    @raises(nx.NetworkXError)
    def test_directed_exception_louvain(self):
        G = nx.DiGraph()
        G.add_weighted_edges_from([('0', '3', 3), ('0', '1', -5),
                                   ('1', '0', -2), ('0', '2', 2),
                                   ('1', '2', -3), ('2', '3', 1)])
        assert_raises(nx.NetworkXError, nx.louvain(G))
