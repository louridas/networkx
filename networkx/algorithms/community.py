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


def louvain(G):
    """Find communities in graph using the Louvain method.

    Parameters
    ----------
    G : NetworkX graph

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
    if G.size(weight='weight') != 0:
        improve = True
        while improve:
            p = _one_pass(G)
            improve, G = _partition_to_graph(p, G)
    else:
        msg = 'The graph has undefined modularity.'
        raise nx.NetworkXError(msg)


def _one_pass(G):
    """
    Puts the Graph nodes into communities as long as it has effect to the
    modularity.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :return: The graph given separated into communities
    """
    increase = True
    p = nx.Graph()
    inner = []
    tot = []
    for u in nx.nodes_iter(G):
        p[u] = u
        _init(G, u, inner, tot)
    while increase:
        increase = False
        for u in nx.nodes_iter(G):
            for v in nx.neighbors(G, u):
                c_old = p[u]
                _remove(u, c_old, p, G, inner, tot)
                c = p[v]
                c_new = p[v]
            _insert(u, c_new, p, G, inner, tot)
            if c_old is not c_new:
                increase = True
    return G


def _partition_to_graph(G, p):
    """
    Creates a Graph that its nodes represent communities and its edges
    represent the connections between these communities.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param p: A partition of the Graph
    :return: If a further Graph partition is possible and a new Graph that has
        a node for each community
    """
    is_possible = True
    return is_possible, G


def _init(G, u, inner, tot):
    """
    Updates the inner and tot parameters.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :param inner: Sum of all the weights of the links inside the community that
        node u is moving into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: Updated inner and tot parameters
    """
    if nx.is_weighted(G, (u, u)):
        inner[u] = G.edge[u][u]['weight']
    else:
        inner[u] = 0
    tot[u] = G.degree(u)


def _remove(G, u, c, p, inner, tot):
    """
    Calculates the inner and tot parameters after extracting the node u from
    the community c.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :param c: The community the node u was before
    :param p: The current partition of the nodes
    :param inner: Sum of all the weights of the links inside the community that
        the node u was moved into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: foo
    """
    if nx.is_weighted(c, (u, u)):
        inner[c] -= _k_in(u, c) + c.edge[u][u]['weight']
    else:
        inner[c] -= _k_in(u, c)
    tot[c] -= G.degree(u)
    p[u] = []


def _insert(G, u, c, p, inner, tot):
    """
    Calculates the inner and tot parameters after inserting the node u into
    the community c.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :param c: The community the node u is moving into
    :param p: The current partition of the nodes
    :param inner: Sum of all the weights of the links inside the community that
        the node u is moving into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: foo
    """
    if nx.is_weighted(c, (u, u)):
        inner[c] += _k_in(u, c) + c.edge[u][u]['weight']
    else:
        inner[c] += _k_in(u, c)
    tot[c] += G.degree(u)
    p[u] = c


def _gain(G, u, c, tot):
    """
    Calculates the change in the modularity.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :param c: The community the node u is going to move into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: The change in modularity
    """
    return float(_k_in(u, c)) / G.size(weight='weight') - float(tot(c) * G.degree(u)) / 2 * pow(G.size(weight='weight'), 2)


def _k_in(u, c):
    """
    Calculates the sum of the weights of the links between node u and other
    nodes in the community.

    It is part of Louvain's algorithm.

    :param u: A node of the Graph
    :param c: The community the node u is moving into
    :return: Sum of the weights of the links between node u and other nodes
        in the community
    """
    community_weight = 0
    for v in nx.neighbors(c, u):
        if nx.is_weighted(c, (u, v)):
            community_weight += c.edge[u][v]['weight']
    return community_weight
