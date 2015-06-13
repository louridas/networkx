# -*- coding: utf-8 -*-
"""Algorithms to detect communities in a Graph."""
import networkx as nx
import networkx.algorithms.isomorphism as iso
__author__ = """\n""".join(['Konstantinos Karakatsanis '
                            '<dinoskarakas@gmail.com>'])
#    Copyright (C) 2004-2015 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
__all__ = ['louvain']


def louvain(G, weight='weight'):
    """Find communities in graph using the Louvain method.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    A dictionary which keys are the nodes of the Graph and its values are the
    communities that they belong

    Notes
    -----
    The Louvain Method for community detection is a method to extract
    communities from large networks created by Vincent Blondel. The method is
    a greedy optimization method that appears to run in time O(n log n).
    """
    H = G.copy()
    tree = []

    if G.number_of_edges() == 0:
        msg = 'The graph has undefined modularity.'
        raise nx.NetworkXError(msg)

    if H.is_directed():
        msg = 'The graph has undefined modularity.'
        raise nx.NetworkXError(msg)

    improve = True
    while improve:
        p = _one_pass(H, weight)
        tree.append(p)
        improve, H = _partition_to_graph(H, p, weight)
    return tree


def _one_pass(G, weight):
    """
    Puts the Graph nodes into communities as long as it has effect to the
    modularity.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :return: The graph given separated into communities
    """
    increase = True
    m = G.size(weight=weight)

    # Initialize p, inner and tot.
    p = {u: u for u in G}
    inner = {u: 0 if not G.has_edge(u, u)
             else G.get_edge_data(u, u).get(weight, 1) for u in G}
    tot = {u: G.degree(u, weight=weight) for u in G}

    # Auxiliary functions.
    k_in = lambda i, com: sum(G.get_edge_data(i, j).get(weight, 1)
                              for j in G.neighbors(i) if p[j] == com)
    gain = lambda n, com: float(k_in(n, com)) / (2 * m) - float(
                tot[com] * G.degree(n, weight=weight)) / (2 * pow(m, 2))

    while increase:
        increase = False
        for u in G:
            c_old = p[u]
            inner[u] -= k_in(u, c_old) if not G.has_edge(u, u) \
                else (k_in(u, c_old)) + G.get_edge_data(u, u).get(weight, 1)
            tot[u] -= G.degree(u, weight=weight)
            max_gain = 0
            best = None
            for v in G.neighbors(u):
                if u != v:
                    c = p[v]
                    diff = gain(u, c)
                    if diff > max_gain:
                        max_gain = diff
                        best = v
            c_new = c_old if best is None else p[best]
            inner[u] += k_in(u, c_new) if not G.has_edge(u, u) \
                else (k_in(u, c_new)) + G.get_edge_data(u, u).get(weight, 1)
            tot[u] += G.degree(u, weight=weight)
            p[u] = c_new
            if c_old != c_new:
                increase = True
    return p


def _partition_to_graph(G, p, weight):
    """
    Creates a Graph that its nodes represent communities and its edges
    represent the connections between these communities.

    It is part of Louvain's algorithm.

    :param G: NetworkX graph
    :param p: A partition of the Graph
    :return: If a further Graph partition is possible and a new Graph that has
        a node for each community
    """
    G2 = nx.Graph()
    G2.add_nodes_from(p.values())
    for u, v, data in G.edges(data=True):
        u_v_weight = data.get(weight, 1)
        edge_weight = u_v_weight if not G2.has_edge(p[u], p[v])\
            else G2.get_edge_data(p[u], p[v]).get(weight, 1) + u_v_weight
        if G2.has_edge(p[u], p[v]):
            G2.edge[p[u]][p[v]][weight] = edge_weight
        else:
            G2.add_edge(p[u], p[v], weight=edge_weight)
    em = iso.numerical_edge_match(weight, 1)
    return (not nx.is_isomorphic(G, G2, edge_match=em)), G2
