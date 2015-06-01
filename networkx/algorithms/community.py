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
            improve, g = _partition_to_graph(p, g)
    else:
        msg = 'The graph has undefined modularity.'
        raise nx.NetworkXError(msg)


def _one_pass(g):
    """
    Puts the Graph nodes into communities as long as it has effect to the
    modularity.

    It is part of Louvain's algorithm.

    :param g: NetworkX graph
    :return: The graph given separated into communities
    """
    increase = True
    p = nx.Graph()
    inner = []
    tot = []
    for i in nx.nodes_iter(g):
        p[i] = i
        _init(g, i, inner, tot)
    while increase:
        increase = False
        for i in nx.nodes_iter(g):
            for j in nx.neighbors(g, i):
                c_old = p[i]
                _remove(i, c_old, p, g, inner, tot)
                c = p[j]
                c_new = p[j]
            _insert(i, c_new, p, g, inner, tot)
            if c_old is not c_new:
                increase = True
    return g


def _partition_to_graph(g, p):
    """
    Creates a Graph that its nodes represent communities and its edges
    represent the connections between these communities.

    It is part of Louvain's algorithm.

    :param g: NetworkX graph
    :param p: A partition of the Graph
    :return: If a further Graph partition is possible and a new Graph that has
        a node for each community
    """
    is_possible = True
    return is_possible, g


def _init(g, i, inner, tot):
    """
    Updates the inner and tot parameters.

    It is part of Louvain's algorithm.

    :param g: NetworkX graph
    :param i: A node of the Graph
    :param inner: Sum of all the weights of the links inside the community that
        node i is moving into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: Updated inner and tot parameters
    """
    if nx.is_weighted(g, (i, i)):
        inner[i] = g.edge[i][i]['weight']
    else:
        inner[i] = 0
    tot[i] = g.degree(i)


def _remove(g, i, c, p, inner, tot):
    """
    Calculates the inner and tot parameters after extracting the node i from
    the community c.

    It is part of Louvain's algorithm.

    :param g: NetworkX graph
    :param i: A node of the Graph
    :param c: The community the node i was before
    :param p: The current partition of the nodes
    :param inner: Sum of all the weights of the links inside the community that
        the node i was moved into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: foo
    """
    if nx.is_weighted(c, (i, i)):
        inner[c] -= _k_in(i, c) + c.edge[i][i]['weight']
    else:
        inner[c] -= _k_in(i, c)
    tot[c] -= g.degree(i)
    p[i] = []


def _insert(g, i, c, p, inner, tot):
    """
    Calculates the inner and tot parameters after inserting the node i into
    the community c.

    It is part of Louvain's algorithm.

    :param g: NetworkX graph
    :param i: A node of the Graph
    :param c: The community the node i is moving into
    :param p: The current partition of the nodes
    :param inner: Sum of all the weights of the links inside the community that
        the node i is moving into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: foo
    """
    if nx.is_weighted(c, (i, i)):
        inner[c] += _k_in(i, c) + c.edge[i][i]['weight']
    else:
        inner[c] += _k_in(i, c)
    tot[c] += g.degree(i)
    p[i] = c


def _gain(g, i, c, tot):
    """
    Calculates the change in the modularity.

    It is part of Louvain's algorithm.

    :param g: NetworkX graph
    :param i: A node of the Graph
    :param c: The community the node i is going to move into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: The change in modularity
    """
    return float(_k_in(i, c)) / g.size(weight='weight') - float(tot(c) * g.degree(i)) / 2 * pow(g.size(weight='weight'), 2)


def _k_in(i, c):
    """
    Calculates the sum of the weights of the links between node i and other
    nodes in the community.

    It is part of Louvain's algorithm.

    :param i: A node of the Graph
    :param c: The community the node i is moving into
    :return: Sum of the weights of the links between node i and other nodes
        in the community
    """
    community_weight = 0
    for j in nx.neighbors(c, i):
        if nx.is_weighted(c, (i, j)):
            community_weight += c.edge[i][j]['weight']
    return community_weight
