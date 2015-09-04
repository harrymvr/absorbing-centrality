from __future__ import division

from itertools import chain

from numpy.random import RandomState
from numpy import arange
from networkx.generators import karate_club_graph, grid_2d_graph
from networkx import is_isomorphic
from absorbing_centrality import keep_largest_component
from absorbing_centrality import is_canonical
from absorbing_centrality import canonical_relabel_nodes
from absorbing_centrality import update_fundamental_matrix
from absorbing_centrality import compute_transition_matrix, compute_fundamental_matrix
from absorbing_centrality import add_supernode
from absorbing_centrality import absorbing_centrality
from absorbing_centrality import SUPER_NODE


def test_label_is_permutation():
    karate_club = karate_club_graph()
    G = canonical_relabel_nodes(karate_club)
    assert set(G.nodes()) == set(range(1, len(G) + 1))
    assert set(G.graph['canonical_map']) ==  set(karate_club.nodes())


def test_is_canonical():
    karate_club = karate_club_graph()
    G = canonical_relabel_nodes(karate_club)
    assert is_canonical(G)
    assert not is_canonical(karate_club)


def test_order_is_consistent():
    karate_club = karate_club_graph()
    G = canonical_relabel_nodes(karate_club)
    for i, v in enumerate(G):
        assert i == v - 1, 'Order of nodes changed after canonicalization.'


def test_largest_component():
    karate_club = karate_club_graph()
    new_club = karate_club.copy()
    # add a triangle, disconnected from the rest of the graph
    new_club.add_edges_from([(100, 101), (101, 102), (102, 100)])
    largest_component = keep_largest_component(new_club)
    assert is_isomorphic(largest_component, karate_club)


def test_update_fundamental_matrix():
    prng = RandomState(20150101)
    P = compute_transition_matrix(karate_club_graph())
    n = P.shape[0]
    order = arange(P.shape[0])
    previous_index = prng.choice(order, 1)
    previous_node = order[previous_index]
    non_absorbing_nodes = chain(range(previous_index),
                                range(previous_index + 1, n))
    non_absorbing_nodes = list(non_absorbing_nodes)
    order = order[non_absorbing_nodes]
    F = compute_fundamental_matrix(P[non_absorbing_nodes, :]
                                   [:, non_absorbing_nodes])
    absorbing_nodes = [previous_node]
    P_updated = P.copy()
    F_updated = F
    while P_updated.shape[0] >= 3:
        next_node = order[prng.choice(len(order), 1)]
        (P_updated, F_updated, order, previous_index) = \
            update_fundamental_matrix(P_updated, F_updated, next=next_node,
                                      previous=previous_node,
                                      previous_index=previous_index,
                                      node_order=order)
        previous_node = next_node
        absorbing_nodes.append(next_node)
        non_absorbing_nodes = [x for x in range(n) if x not in absorbing_nodes]
        F_slow = compute_fundamental_matrix(P[non_absorbing_nodes, :]
                                            [:, non_absorbing_nodes])
        error_at_step = sum(sum(F_updated - F_slow).T)[0, 0]
        assert abs(error_at_step) < 1e-8, "Error is more than 1e-8."


def test_add_supernode():
    karate_club = karate_club_graph()
    n = len(karate_club)
    m = len(karate_club.edges())
    karate_club = add_supernode(karate_club, query=[10, 20, 30])
    assert len(karate_club) == n + 1, 'Supernode was not added.'
    # for v in [10, 20, 30]:
        # assert(n + 1,
               # karate_club.graph['canonical_map'][v] in karate_club.edges(),
        # 'Edge from supernode to query not added.')
    assert len(karate_club.edges()) == 2 * m + 3, \
        'Number of edges after adding supernode is incorrect.'


def test_supernode_last():
    karate_club = karate_club_graph()
    karate_club_super = add_supernode(karate_club, query=[])
    for i, v in enumerate(karate_club_super):
        if i == len(karate_club_super) - 1:
            assert karate_club_super.graph['label_map'][v] == SUPER_NODE


def test_transition_matrix_with_supernode():
    karate_club = karate_club_graph()
    karate_club = canonical_relabel_nodes(keep_largest_component(karate_club))
    karate_club_super = add_supernode(karate_club, query=[10, 20, 30])
    P = compute_transition_matrix(karate_club_super)
    for v in [10, 20, 30]:
        assert P[-1, v - 1] == 1 / 3


def test_absorbing_centrality():
    dimension = 4
    grid = grid_2d_graph(dimension, dimension, periodic=False).to_undirected()
    solutions = {
        (0, 0): 34.7857142857,
        (0, 1): 22.9285714286,
        (0, 2): 22.9285714286,
        (0, 3): 34.7857142857,
        (1, 0): 22.9285714286,
        (1, 1): 14.0714285714,
        (1, 2): 14.0714285714,
        (1, 3): 22.9285714286,
        (2, 0): 22.9285714286,
        (2, 1): 14.0714285714,
        (2, 2): 14.0714285714,
        (2, 3): 22.9285714286,
        (3, 0): 34.7857142857,
        (3, 1): 22.9285714286,
        (3, 2): 22.9285714286,
        (3, 3): 34.7857142857
    }
    for i in range(dimension):
        for j in range(dimension):
            assert abs(absorbing_centrality(grid, [(i, j)], epsilon=1e-12) -
                       solutions[(i, j)]) < 1e-5
