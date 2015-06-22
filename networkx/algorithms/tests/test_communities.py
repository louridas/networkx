#!/usr/bin/env python
from nose.tools import *
import networkx as nx


class TestCommunity:

    def test_louvain_undirected(self):
        G = nx.Graph()
        G.add_edges_from([(0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 4),
                          (1, 7), (2, 4), (2, 5), (2, 6), (3, 7), (4, 10),
                          (5, 7), (5, 11), (6, 7), (6, 11), (8, 9), (8, 10),
                          (8, 11), (8, 14), (8, 15), (9, 12), (9, 14), (10, 11),
                          (10, 12), (10, 13), (10, 14), (11, 13)])
        tree = nx.louvain(G)
        assert_equal(len(tree), 3)
        t1 = tree[0]
        assert_equal(t1.nodes(data=True),
                     [(4, {'members': set([0, 1, 2, 4, 5])}),
                      (12, {'members': set([8, 9, 10, 12, 14, 15])}),
                      (13, {'members': set([11, 13])}),
                      (7, {'members': set([3, 6, 7])})])
        assert_equal(t1.edges(data=True),
                     [(4, 12, {'weight': 1}),
                      (4, 4, {'weight': 14}),
                      (4, 13, {'weight': 1}),
                      (4, 7, {'weight': 4}),
                      (12, 12, {'weight': 16}),
                      (12, 13, {'weight': 3}),
                      (13, 13, {'weight': 2}),
                      (13, 7, {'weight': 1}),
                      (7, 7, {'weight': 4})])
        t2 = tree[1]
        assert_equal(t2.nodes(data=True),
                     [(2, {'members': set([12, 13])}),
                      (3, {'members': set([4, 7])})])
        assert_equal(t2.edges(data=True),
                     [(2, 2, {'weight': 24}),
                      (2, 3, {'weight': 3}),
                      (3, 3, {'weight': 26})])
        t3 = tree[2]
        assert_equal(t3.nodes(data=True),
                     [(1, {'members': set([2, 3])})])
        assert_equal(t3.edges(data=True),
                     [(1, 1, {'weight': 56})])

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

