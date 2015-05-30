# -*- coding: utf-8 -*-
"""Algorithms to detect communities in a Graph."""
import networkx as nx
__author__ = """\n""".join(['Konstantinos Karakatsanis <dinoskarakas@gmail.com>',
                            'Thodoris Sotiropoulos <theosotr@windowslive.com>'])
#    Copyright (C) 2004-2015 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
__all__ = ['louvain']


def louvain(g):
    """Find communities in graph using the Louvain method.

    Parameters
    ----------
    g : NetworkX graph

    Returns
    -------
    foo

    Examples
    --------
    >>> G = nx.path_graph(10)
    >>> comp = louvain(G)
    >>> comp[0]
    foo

    Notes
    -----
    The Louvain Method for community detection is a method to extract
    communities from large networks created by Vincent Blondel. The method is
    a greedy optimization method that appears to run in time O(n log n).
    """
    if g.size(weight='weight') != 0:
        improve = True
        while improve:
            p = _one_pass(g)
            g = _partition_to_graph(p, g)
            improve = False
    else:
        msg = 'The graph has undefined modularity.'
        raise nx.NetworkXError(msg)


def _one_pass(g):
    increase = True
    p = nx.Graph()
    for i in nx.nodes_iter(g):
        p[i] = i
        _init(g, i)
    while increase:
        increase = False
        for i in nx.nodes_iter(g):
            for j in nx.neighbors(g, i):
                c_old = p[i]
                _remove(i, c_old, p)
                c = p[j]
                c_new = p[j]
            _insert(i, c_new, p)
            if c_old is not c_new:
                increase = True
    return g


def _partition_to_graph(p, g):
    return p, g


def _init(g, i):
    if nx.is_weighted(g, (i, i)):
        init[i] = g.edge[i][i]['weight']
    else:
        init[i] = 0
    tot[i] = g.degree(i)


def _remove(i, c, p):
    if nx.is_weighted(c, (i, i)):
        init[c] -= 2 * _k_in(i, c) + c.edge[i][i]['weight']
    else:
        init[c] -= 2 * _k_in(i, c)
    tot[c] -= g.degree(i)
    p[i] = []


def _insert(i, c, p):
    if nx.is_weighted(c, (i, i)):
        init[c] += 2 * _k_in(i, c) + c.edge[i][i]['weight']
    else:
        init[c] += 2 * _k_in(i, c)
    tot[c] -= g.degree(i)
    tot[c] += g.degree(i)
    p[i] = c


def _gain(i, c, g):
    return float(_k_in(i, c)) / g.size(weight='weight') - float(_tot(c) * g.degree(i)) / 2 * pow(g.size(weight='weight'), 2)


def _k_in(i, c):
    community_weight = 0
    for j in nx.neighbors(c, i):
        if nx.is_weighted(c, (i, j)):
            community_weight += c.edge[i][j]['weight']
    return community_weight
