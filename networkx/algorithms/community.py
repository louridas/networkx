# -*- coding: utf-8 -*-
"""Algorithms to detect communities in a Graph."""
import networkx as nx
import networkx.algorithms.isomorphism as iso
__author__ = """\n""".join(['Konstantinos Karakatsanis <dinoskarakas@gmail.com>'])
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
    if G.size(weight='weight') == 0:
        msg = 'The graph has undefined modularity.'
        raise nx.NetworkXError(msg)
    elif G.is_directed():
        msg = 'The graph has undefined modularity.'
        raise nx.NetworkXError(msg)
    improve = True
    while improve:
        p = _one_pass(G)
        improve, G = _partition_to_graph(G, p)


def _one_pass(G):
    """
    Puts the Graph nodes into communities as long as it has effect to the
    modularity.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :return: The graph given separated into communities
    """
    increase = True
    p = {}
    inner = {}
    tot = {}
    for u in G:
        p[u] = u
        inner[u], tot[u] = _init(G, u)
    while increase:
        increase = False
        for u in G:
            c_old = p[u]
            p, inner[u], tot[u] = _remove(G, u, c_old, p, inner[u], tot[u])
            if G.neighbors(u) != 0:
                max_gain = 0
                for v in nx.neighbors(G, u):
                    if u != v:
                        c = p[v]
                        if _gain(G, u, c, p, tot[c]) > max_gain:
                            max_gain = _gain(G, u, c, p, tot[c])
                            best = v
            try:
                c_new = p[best]
            except NameError:
                c_new = c_old
            p, inner[u], tot[u] = _insert(G, u, c_new, p, inner[u], tot[u])
            if c_old is not c_new:
                increase = True
    return p


def _partition_to_graph(G, p):
    """
    Creates a Graph that its nodes represent communities and its edges
    represent the connections between these communities.

    It is part of Louvain's algorithm.

    After https://goo.gl/BB78Mv lines 304-312

    :param G: NetworkX graph
    :param p: A partition of the Graph
    :return: If a further Graph partition is possible and a new Graph that has
        a node for each community
    """
    G2 = nx.Graph()
    G2.add_nodes_from(p.values())
    for u, v, data in G.edges(data=True):
        u_v_weight = data.get("weight", 1)
        weight = G2.get_edge_data(p[u], p[v], {'weight': 0}).get('weight', 1) + u_v_weight
        G2.add_edge(p[u], p[v], weight=weight)
    em = iso.numerical_edge_match('weight', 1)
    return (not nx.is_isomorphic(G, G2, edge_match=em)), G2


def _init(G, u):
    """
    Updates the inner and tot parameters.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :return: Updated inner and tot parameters
    """
    inner = G.get_edge_data(u, u, {'weight': 0}).get('weight', 1)
    tot = G.degree(u)
    return inner, tot


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
    :return: Updated p, inner and tot parameters
    """
    inner -= _k_in(G, u, c, p) + G.get_edge_data(u, u, {'weight': 0}).get('weight', 1)
    tot -= G.degree(u)
    p[u] = []
    return p, inner, tot


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
    :return: Updated p, inner and tot parameters
    """
    inner += _k_in(G, u, c, p) + G.get_edge_data(u, u, {'weight': 0}).get('weight', 1)
    tot += G.degree(u)
    p[u] = c
    return p, inner, tot


def _gain(G, u, c, p, tot):
    """
    Calculates the change in the modularity.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :param c: The community the node u is going to move into
    :param tot: Sum of all the weights of the links to nodes in the community
    :return: The change in modularity
    """
    m = G.size(weight='weight')
    return float(_k_in(G, u, c, p)) / (2 * m) - float(tot * G.degree(u)) / (2 * pow(m, 2))


def _k_in(G, u, c, p):
    """
    Calculates the sum of the weights of the links between node u and other
    nodes in the community.

    It is part of Louvain's algorithm.

    :param u: A node of the Graph
    :param c: The community the node u is moving into
    :return: Sum of the weights of the links between node u and other nodes
        in the community
    """
    return sum(G.get_edge_data(u, v, {'weight': 0}).get('weight', 1) for v in G if p[v] == c)
