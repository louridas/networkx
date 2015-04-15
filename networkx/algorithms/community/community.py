#-*- coding: utf-8 -*-
import networkx as nx

__all__ = ['girvan_newman']


def girvan_newman(G, weight=None):
    """Find communities in graph using Girvan-Newman method.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default=None)
       Edge data key corresponding to the edge weight.

    Returns
    -------
    List of tuples which contains the clusters of nodes.

    Examples
    --------
    >>> G = nx.path_graph(10)
    >>> comp = girvan_newman(G)
    >>> print(comp[0])
    ([0, 1, 2, 3, 4], [8, 9, 5, 6, 7])

    Notes
    -----
    The Girvan–Newman algorithm detects communities by progressively removing
    edges from the original graph. Algorithm removes edge with the highest
    betweenness centrality at each step. As the graph breaks down into pieces,
    the tightly knit community structure is exposed and result can be depicted
    as a dendogram.
    """
    g = G.copy().to_undirected()
    components = []
    while g.number_of_edges() > 0:
        _remove_max_edge(g, weight)
        components.append(_communities(g),)
    return components


def _remove_max_edge(G, weight=None):
    """
    Removes edge with the highest value on betweenness centrality.

    Repeat this step until more connected components than the connected
    components of the original graph are detected.

    It is part of Girvan-Newman algorithm.

    :param G: NetworkX graph
    :param weight: string, optional (default=None) Edge data key corresponding
    to the edge weight.
    """
    number_components = nx.number_connected_components(G)
    while nx.number_connected_components(G) <= number_components:
        betweenness = nx.edge_betweenness_centrality(G, weight=weight)
        max_value = max(betweenness.values())
        for edge in G.edges():
            if betweenness[edge] == max_value:
                G.remove_edge(*edge)


def _communities(G):
    """
    Stores nodes of every connected component of graph in a tuple.

    It is part of Girvan-Newman algorithm.

    :param G: NetworkX graph
    :return: Tuple which contains the list of nodes of every connected component
    detected in graph.
    """
    communities_comp = ()
    for subgraph in nx.connected_component_subgraphs(G):
        communities_comp += (subgraph.nodes(),)
    return communities_comp
