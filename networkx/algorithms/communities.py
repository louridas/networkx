# -*- coding: utf-8 -*-
"""Algorithms to detect communities in a Graph."""
import networkx as nx
from collections import defaultdict
__author__ = """\n""".join(['Konstantinos Karakatsanis '
                            '<dinoskarakas@gmail.com>',
                            'Panos Louridas <louridas@gmail.com>'])
#    Copyright (C) 2015 by
#    Konstantinos Karakatsanis <dinoskarakas@gmail.com>
#    Panos Louridas <louridas@gmail.com>
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

    improved = True
    while improved:
        p, c, improvements = _one_pass(H)
        H = _partition_to_graph(H, p, c)
        improved = improvements > 0 and len(c) > 1
        tree.append(H)
    return tree

def _one_pass(G):
    """
    A single pass of the Louvain algorithm. Moves nodes from community
    to community as long as this increases the graph's modularity.

    Part of the Louvain algorithm.

    :param G: NetworkX graph
    :return: a tuple (p, c, improvements); p is a dictionary whose keys
        are the nodes of the graph and the values are the communities
        they belong; c is a dictionary whose keys are the communities and
        the values are the nodes they contain; improvements is the number
        of times the modularity was improved by moving a node to 
        another community
    """
    i = 0 # community index, starting from zero
    # partition dict: each entry has as key a node index and as value
    # the index of the community in which the node belongs
    p = {}
    # community dict: each entry has as key a community index and as
    # value a set containing the indices of the nodes that belong in the
    # community
    c = {}
    # Sum of all the weights of the links inside each community
    inner = {} # the number of internal links of a community
    # Sum of all the weights of the links to nodes in each community
    tot = {} 
    # The denominator used to calculate the gain in modularity;
    # we precompute it and cache it to avoid re-calculations in loop
    # calls.
    denom = 2 * G.size(weight='weight')
    improvements = 0 # number of improvements
    for u in G:
        p[u] = i
        c[i] = { u }
        inner[i], tot[i] = _init(G, u)
        i += 1
    increase = True        
    while increase:
        increase = False
        for u in G:
            c_old = p[u]
            _remove(G, u, c_old, c, p, inner, tot)
            if G.neighbors(u) != 0:
                max_gain = 0
                c_new = None
                for v in nx.neighbors(G, u):
                    if u != v:
                        c_v = p[v]
                        gain = _gain(G, u, c_v, p, tot, denom)
                        if gain > max_gain:
                            max_gain = gain
                            c_new = c_v
            if c_new is None:
                c_new = c_old
            _insert(G, u, c_new, c, p, inner, tot)
            if c_new != c_old:
                increase = True
                improvements += 1
    c = { key:value for (key, value) in c.iteritems() if len(value) > 0 }
    return p, c, improvements

def _partition_to_graph(G, p, c):
    """Creates a graph whose nodes represent communities and its edges
    represent the connections between these communities.

    Part of the Louvain algorithm.

    :param G: NetworkX graph
    :param p: Partinion dictionary; each entry has as key a node index
        and as value the index of the community in which the node belongs    
    :param c: Community dictionary; each entry has as key a community index
        and as value a set containing the indices of the nodes that belong
        in the community

    :return: a new Graph that has a node for each community

    """
    G2 = nx.Graph()
    for k in c:
        G2.add_node(k, members=c[k])
    
    for u, v, data in G.edges(data=True):
        u_v_weight = data.get('weight', 1)
        if p[u] == p[v] and u != v:
            u_v_weight *= 2
        if G2.has_edge(p[u], p[v]):
            G2.edge[p[u]][p[v]]['weight'] += u_v_weight
        else:
            G2.add_edge(p[u], p[v], weight=u_v_weight)
    return G2

def _init(G, u):
    """
    Returns a tuple (inner, tot) where inner is the sum of weights of
    links from u to u, i.e., loops, and tot is the sum of weights of all
    links adjacent to u, i.e., its degree.

    Part of the Louvain algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :return: Updated inner and tot parameters
    """
    inner = G.get_edge_data(u, u, {'weight': 0}).get('weight', 1)
    tot = G.degree(u)
    return inner, tot


def _remove(G, u, c_old, c, p, inner, tot):
    """
    Adjusts the p, inner, and tot parameters after extracting the node u from
    the community c.

     Part of the Louvain algorithm.

    :param G: NetworkX graph
    :param u: A node of the graph
    :param c_old: The community that contains node u
    :param c: The current communities of the nodes
    :param p: The current partition of the nodes; when the function returns p
        has been removed from its partition, so p[u] is None
    :param inner: Sum of all the weights of the links inside each community;
        when the function returns the weights of the internal edges
        adjacent to u have been removed
    :param tot: Sum of all the weights of the links to nodes in each community;
        when the function returns all the weights of edges adjacent to u
        have been removed
    """
    inner[c_old] -= _k_in(G, u, c, p) + \
                    G.get_edge_data(u, u, {'weight': 0}).get('weight', 1)
    tot[c_old] -= G.degree(u)
    p[u] = None
    c[c_old].remove(u)

def _insert(G, u, c_new, c, p, inner, tot):
    """
    Calculates the inner and tot parameters after inserting the node u into
    the community c.

    Part of the Louvain algorithm.

    :param G: NetworkX graph
    :param u: A node of the Graph
    :param c_new: The community the node u is moving into
    :param c: The current communities of the nodes
    :param p: The current partition of the nodes; when the function returns
        we have p[u] = c_new
    :param inner: Sum of all the weights of the links inside each community;
        when the function returns the weights of the internal edges
        adjacent to u have been removed
    :param tot: Sum of all the weights of the links to nodes in the community;
        when the function returns all the weights of edges adjacent to u
        have been removed       
    :return: Updated p, inner and tot parameters
    """
    inner[c_new] += _k_in(G, u, c, p) + \
                    G.get_edge_data(u, u, {'weight': 0}).get('weight', 1)
    tot[c_new] += G.degree(u)
    c[c_new].add(u)
    p[u] = c_new

def _gain(G, u, c_target, p, tot, denom):
    """
    Calculates the gain in modularity.

    Part of the Louvain algorithm.

    :param G: NetworkX graph
    :param u: A node of the graph
    :param c_target: The community the node u is going to move into
    :param tot: Sum of all the weights of the links to nodes in each community
    :param tot: Sum of all the weights of the links to nodes in the communit    
    :param denom: Twice the sum of all the weights of all links.
    :return: The change in modularity
    """
    return float(_k_in(G, u, c_target, p)) - \
           float(tot[c_target] * G.degree(u)) / (denom)


def _k_in(G, u, c_target, p):
    """
    Calculates the sum of the weights of the links between node u and other
    nodes in the community.

    Part of the Louvain algorithm.

    :param u: A node of the graph
    :param c_target: The community the node u is moving into
    :return: Sum of the weights of the links between node u and other nodes
        in the community
    """
    return sum(G.get_edge_data(u, v, {'weight': 0}).get('weight', 1)
               for v in G if v != u and p[v] == c_target)

