"""
    graph_helpers
    ~~~~~~~~~~

    Provides helper functions that are used for:
        * graph preprocessing
        * matrix genereration (fundamental, transition)
        * Sherman-Morrison formula
        * absorbing random-walk centrality computation


    :copyright: 2015 Charalampos Mavroforakis, <cmav@bu.edu> and contributors
    :license: ISC

"""
from __future__ import division

from itertools import chain

from networkx import (DiGraph, adjacency_matrix, connected_component_subgraphs,
                      relabel_nodes)
from numpy import arange, argwhere, array, eye, matrix, zeros
from scipy.sparse import eye as speye, csc_matrix, diags
from scipy.sparse.linalg import inv as spinv, spilu

from .exceptions import CanonicalizationError

SUPER_NODE = '_super_'

try:
    dict.iteritems
except AttributeError:
    # Python 3
    def itervalues(d):
        return iter(d.values())

    def iteritems(d):
        return iter(d.items())
else:
    # Python 2
    def itervalues(d):
        return d.itervalues()

    def iteritems(d):
        return d.iteritems()


def keep_largest_component(G):
    """Keeps the largest connected component of the graph and removes all the
    self loops.

    Parameters
    ----------
    G : NetworkX graph
        The input graph.

    Returns
    -------
    NetworkX graph
        The largest connected subgraph of G. All self-loops have been removed.

    """
    G_largest = sorted(connected_component_subgraphs(G),
                       key=len, reverse=True)[0]
    G_largest.remove_edges_from(G_largest.selfloop_edges())
    return G_largest


def is_canonical(G):
    """Tests if the graph has been canonicalized.

    Parameters
    ----------
    G : NetworkX graph


    Returns
    -------
    bool
    Returns True, if the graph has been canonicalized.

    """
    if 'label_map' in G.graph and 'canonical_map' in G.graph:
        return True
    else:
        return False


def canonical_relabel_nodes(G):
    """Relabels the nodes in the graph, such that the new names belong in
    the set [1,n]. The labeling information is stored in the dictionaries
    `G.graph['canonical_map']` and `G.graph['label_map']`. These provide a way
    to map original to canonical node names and vice-versa, respectively.

    Parameters
    ----------
    G : NetworkX graph

    Returns
    -------
    G_prime : NetworkX graph
        The relabeled graph. It includes two attributes:

        (i) `G_prime.graph['canonical_map']` : dict
                Holds the mapping between the original names of the nodes and
                the new, canonical, names (original -> new).

        (ii) `G_prime.graph['label_map']` : dict
                Holds the mapping between the new, canonical, names and the
                original names of the nodes (new -> original).


    Note: The relabeling of a particular node might not be consistent
    across two consecutive runs. Also, the relabeling happens on a copy, so the
    original graph will be untouched.

    """
    n = G.number_of_nodes()
    mapping = dict(zip(G.nodes(), range(1, n + 1)))
    G_prime = relabel_nodes(G, mapping, copy=True)
    G_prime.graph['canonical_map'] = mapping.copy()
    G_prime.graph['label_map'] = dict((v, k) for k, v in iteritems(mapping))
    return G_prime


def compute_fundamental_matrix(P, fast=True, drop_tol=1e-5, fill_factor=1000):
    """Computes the fundamental matrix for an absorbing random walk.

    Parameters
    ----------
    P : scipy.sparse matrix
        The transition probability matrix of the absorbing random walk. To
        construct this matrix, you start from the original transition matrix
        and delete the rows that correspond to the absorbing nodes.

    fast : bool, optional
    If True (default), use the iterative SuperLU solver from
    scipy.sparse.linalg.

    drop_tol : float, optional
        If `fast` is True, the `drop_tol` parameter of the SuperLU solver is
        set to this value (default is 1e-5).

    fill_factor: int, optional
        If `If `fast` is True, the `fill_factor` parameter of the SuperLU
        solver is set to this value (default is 1000).

    Returns
    -------
    F : scipy.sparse matrix
        The fundamental matrix of the random walk. Element (i,j) holds the
        expected number of times the random walk will be in state j before
        absorption, when it starts from state i. For more information, check
        [1]_.

    References
    ----------
    .. [1] Doyle, Peter G., and J. Laurie Snell.
       Random walks and electric networks.
       Carus mathematical monographs 22 (2000).
       https://math.dartmouth.edu/~doyle/docs/walks/walks.pdf

    """
    n = P.shape[0]
    F_inv = speye(n, format='csc') - P.tocsc()
    if fast:
        solver = spilu(F_inv, drop_tol=drop_tol, fill_factor=fill_factor)
        F = matrix(solver.solve(eye(n)))
    else:
        F = spinv(F_inv).todense()
    return F


def compute_personalized_transition_matrix(G, alpha=0.85,
                                           restart_set=[SUPER_NODE]):
    """Returns the transition matrix of the random walk with restarts.

    Parameters
    ----------
    G : graph

    alpha : float, optional
        The probability of the random surfer to continue their walk (default is
        0.85).

    restart_set : list, optional
        The set of nodes to restart from. If not supplied, the restarts lead to
        the supernode (default is [SUPER_NODE]).


    Returns
    -------
    P : scipy.sparse.matrix
        The probability matrix for the random walk with restarts.

    """
    if not has_supernode(G) and SUPER_NODE in restart_set:
        raise CanonicalizationError('Cannot restart the random walks at the '
                                    'supernode')
    canonical_restart_set = [G.graph['canonical_map'][n] for n in restart_set]
    restart_graph = DiGraph()
    restart_edges = [(n, q) for n in G.nodes() for q in canonical_restart_set]
    # add a self loop edge at the supernode, if there is one, to avoid division
    # by zero when computing the transition matrix
    if has_supernode(G):
        # TODO Why is this edge added -- and then removed?
        restart_edges.append((G.graph['canonical_map'][SUPER_NODE],
                              G.graph['canonical_map'][SUPER_NODE]))
    restart_graph.add_edges_from(restart_edges)
    P = compute_transition_matrix(G)
    P_restart = compute_transition_matrix(restart_graph)
    if has_supernode(G):
        # remove the bottom right corner (the added self-loop)
        P_restart[-1, -1] = 0
    P_final = alpha * P + (1 - alpha) * P_restart
    # TODO do the transition probabilities from SUPER_NODE sum up to 1?
    # Does it matter?
    return P_final


def compute_transition_matrix(G):
    r"""Builds the random transition matrix P. The probability of going from
    node `i` to node`j` is equal to:

    .. math::
       P_{i,j} = \frac{1}{\text{degree}(i)}

    Parameters
    ----------
    G : NetworkX graph


    Returns
    -------
    P : scipy.sparse matrix
        The random transition probability matrix.

    """
    # TODO: make this function work for weighted graphs
    adjacency = adjacency_matrix(G)
    degrees = adjacency.sum(axis=1)
    inverse_degrees = 1 / degrees
    P = diags(inverse_degrees.T.tolist(), [0]) * adjacency
    return P


def _fast_update_fundamental_rows(P, F, row=0, row_previous=0):
    """Updates the rows of the fundamental matrix by adding an absorbing node,
    using the Woodbury identity.

    Parameters
    ----------
    P : matrix
        The transition matrix of the graph.

    F : matrix
        The fundamental matrix of the graph after setting `row_previous` to be
        an absorbing node.

    row : int, optional
        The index of the row of `F` that will be set to be absorbing (default is
        0).

    row_previous : int, optional
        The index of the row of `P` that was set to be absorbing when
        computing `F` (default is 0).

    Returns
    -------
    F_updated_row : matrix
        The new fundamental matrix, where row of the absorbing node has been
        updated.

    """
    n_P = P.shape[0]
    n_F = F.shape[0]
    # These are the rows of the matrix P that were used in order to compute F
    previous_non_absorbing = chain(range(row_previous),
                                   range(row_previous + 1, n_P))
    previous_non_absorbing = list(previous_non_absorbing)
    # Re-order the transition matrix P without row and column `row_previous`,
    # such that the row (and column) that is absorbing will be first.
    reordering = [row] + list(chain(range(row),
                                    range(row + 1, n_F)))
    P_non_absorbing = P[previous_non_absorbing, :][:, previous_non_absorbing]
    F_reordered = F.take(reordering, axis=0).take(reordering, axis=1)
    U = csc_matrix((n_F, 1))
    U[0, 0] = 1.0
    V = P_non_absorbing[row, reordering]
    y = 1.0 + V.dot(F_reordered[:, 0])[0, 0]
    VF = V.dot(F_reordered)
    F_updated_row = F_reordered - F_reordered[:, 0].dot(VF) / y
    return F_updated_row


def _fast_update_fundamental_columns(P, F, col=0, col_previous=0):
    """Updates the columns of the fundamental matrix by adding an absorbing
    node, using the Woodbury identity.

    Parameters
    ----------
    P : matrix
        The transition matrix of the graph.

    F : matrix
        The fundamental matrix of the graph after setting `col_previous` to be
        an absorbing node. `F` needs to be the result of the
        _fast_update_fundamental_rows(), otherwise the order of the columns is
        wrong.

    col : int, optional
        The index of the column of `F` that corresponds to the new absorbing
        node (default is 0).

    col_previous : int, optional
        The index of the column of `P` that corresponds to the absorbing node
        when computing `F` (default is 0).

    Returns
    -------
    F_updated_col : matrix
        The new fundamental matrix, in which the column of the absorbing
        node has been updated.

    """
    n_P = P.shape[0]
    n_F = F.shape[0]

    # These are the rows of the matrix P that were used in order to compute FP
    previous_round_indices = chain(range(col_previous),
                                   range(col_previous + 1, n_P))
    previous_round_indices = list(previous_round_indices)
    # Re-order the transition matrix P without row and column `col_previous`,
    # such that the col (and row) that is absorbing will be first.
    current_indices = [col] + list(chain(range(col),
                                         range(col + 1, n_F)))
    P_non_absorbing = P[previous_round_indices, :][:, previous_round_indices]
    # F is the result of the _fast_update_fundamental_rows(), so it is already
    # re-arranged in the way we want it.
    F_reordered = F
    U = P_non_absorbing[current_indices, col].A.T
    y = 1.0 + F_reordered[0, :].dot(U)[0, 0]
    FU = F_reordered.dot(U)
    F_updated_col = F_reordered - FU.dot(F_reordered[0, :]) / y
    return F_updated_col


def update_fundamental_matrix(P, F, next, previous, previous_index=0,
                              node_order=None):
    """Applies Woodbury's formula to update the fundamental matrix in order to
    avoid doing an inversion.

    Parameters
    ----------
    P : matrix
        The transition matrix of the graph, where `previous` is non absorbing.

    F : matrix
        The fundamental matrix of the graph after setting the node `previous`
        as an absorbing node.

    next : int
        The node that will be set as absorbing next. The result of this call
        will result in a fundamental matrix where `next` is an absorbing node.

    previous : int
        The node that was set as absorbing when computing F.

    previous_index : int, optional
        The row/col index of node `previous` in P (default is 0).

    node_order : list, optional
        The nodes that corresponds to the rows/cols of `P`, in order. If not
        supplied, the order is considered to be [0, .. , n_P - 1], where
        n_P is the the number of rows/cols in `P`.


    Returns
    -------
    P_updated : matrix
        The new transition matrix, where `previous` is absorbing.

    F_updated : matrix
        The fundamental matrix after adding `previous` and `next` in the set of
        absorbing nodes.

    node_order_updated : list
        The new order of the non absorbing nodes in the `F_new`.

    next_index : int
        The row/col index of the node `next`, that we just set as absorbing,
        in `P_new`.

    """
    n_P = P.shape[0]
    n_F = F.shape[0]
    if node_order is None:
        node_order = arange(n_P)
        node_order = node_order[list(chain(range(previous),
                                           range(previous + 1, n_P)))]
        previous_index = previous
    if type(node_order) is list:
        node_order = array(node_order)
    # We need to find which row/column corresponds to the node `next`
    next_index = argwhere(array(node_order) == next)[0]
    node_order = node_order[list(chain(range(next_index),
                                       range(next_index + 1, n_F)))]
    row_update = _fast_update_fundamental_rows(P, F, row=next_index,
                                               row_previous=previous_index)
    col_update = _fast_update_fundamental_columns(P, row_update,
                                                  col=next_index,
                                                  col_previous=previous_index)
    P_updated = P[list(chain(range(previous_index),
                             range(previous_index + 1, P.shape[0]))), :] \
                 [:, list(chain(range(previous_index),
                                range(previous_index + 1, P.shape[0])))]
    F_updated = col_update[1:, 1:]
    node_order_updated = node_order.tolist()
    return P_updated, F_updated, node_order_updated, next_index


def add_supernode(G, query=None):
    """Adds a supernode to the graph and connects it with directed edges to the
    query nodes.

    Parameters
    ----------
    G : NetworkX graph
        The graph in which we want to add a supernode.

    query : list, default is None
        The list of nodes that the supernode will be connected to. If `query`
        is None, the supernode will be connected to all the nodes in `G`. These
        new edges will be directed.


    Returns
    -------
    NetworkX graph
        A directed graph with the supernode, and the new edges, added. The
        attributes of the graph, i.e. 'label_map' and 'canonical_map', are also
        updated (or created if the input graph was not canonicalized) to
        reflect the new node.

    """
    if has_supernode(G):
        return G
    n = len(G)
    if query is None or not query:
        query = G.nodes()
    if not is_canonical(G):
        G = canonical_relabel_nodes(keep_largest_component(G))
        n = len(G)
        query = [G.graph['canonical_map'][q] for q in query]
    if not G.is_directed():
        G = G.to_directed()
    G.add_node(n + 1)
    if SUPER_NODE in G.graph['canonical_map']:
        # make sure that SUPER_NODE is the last node returned by nodes()
        nodes = G.nodes()
        if not G.graph['label_map'][nodes[-1]] == SUPER_NODE:
            raise CanonicalizationError('Supernode is not the last node.')
    G.graph['canonical_map'][SUPER_NODE] = n + 1
    G.graph['label_map'][n + 1] = SUPER_NODE
    G.add_edges_from([(n + 1, q) for q in query])
    return G


def has_supernode(G):
    """Checks if there exist a supernode in the graph.

    Parameters
    ----------
    G : NetworkX graph


    Returns
    -------
    has_supernode : bool

    """
    if is_canonical(G) and SUPER_NODE in G.graph['canonical_map']:
        return True
    else:
        return False


def absorbing_centrality_inversion(G, team, query=None, with_restarts=False,
                                   alpha=0.85):
    """Compute the absorbing centrality of a team using a fast inversion with
    SuperLU solver.

    Parameters
    ----------
    G : NetworkX graph
        The graph on which to compute the centrality.

    team : list
        The team of nodes, whose centrality to compute.

    query : list, optional
        The set of query nodes to use for the random walks. If None (default)
        or empty, the query set is equal to the set of all nodes in the graph.

    with_restarts : bool, optional
        If True, restarts the random surfer to the the query set (default is
        False).

    alpha : float, optional
        The probability of the random surfer to continue (default is 0.85).


    Returns
    -------
    score : float
        The absorbing centrality score.


    Note
    ----
    Both `team` and `query` should  use the original node names.

    """
    if not is_canonical(G):
        G = canonical_relabel_nodes(keep_largest_component(G))
    if query is None:
        query = []
    canonical_query = [G.graph['canonical_map'][q] for q in query]
    G = add_supernode(G, canonical_query)
    if with_restarts:
        P = compute_personalized_transition_matrix(G, alpha)
    else:
        P = compute_transition_matrix(G)
    team_ind = [G.graph['canonical_map'][c] for c in team]
    non_absorbing_nodes = [v - 1 for v in G if v not in team_ind]
    P_non_absorbing = P[non_absorbing_nodes, :][:, non_absorbing_nodes]
    F = compute_fundamental_matrix(P_non_absorbing)
    row_sums = F.sum(axis=1)
    score = row_sums[-1].sum()
    return score


def absorbing_centrality(G, team, query=None, P=None, epsilon=1e-5,
                         max_iterations=None, with_restarts=False,
                         alpha=0.85):
    r"""Compute the absorbing centrality of a team. The algorithm works by
    iteratively computing the powers of the non-absorbing submatrix of the
    transition matrix P.

    Parameters
    ----------
    G : NetworkX graph
        The graph on which to compute the centrality.

    team : list
        The team of nodes, whose centrality to compute.

    query : list, optional
        The set of query nodes to use for the random walks. If None (default)
        or empty, the query set is equal to the set of all nodes in the graph.

    P : matrix, optional
        The precomputed transition matrix of the graph (default is None).

    epsilon : float, optional
        The iterative algorithm stops when the error between the centrality
        computed by two successive iterations falls below epsilon (default is
        1e-5).

    max_iterations : int, optional
        The upper limit to the number of iteratios of the algorithm (default
        is None).

    with_restarts : bool, optional
        If True, restarts the random surfer to the the query set (default is
        False).

    alpha : float, optional
        The probability of the random surfer to continue (default is 0.85).


    Returns
    -------
    score : float
        The absorbing centrality score.


    Note
    ----
    Both `team` and `query` should  use the original node names.

    """
    if not is_canonical(G):
        G = canonical_relabel_nodes(keep_largest_component(G))
    if query is None:
        query = []
    canonical_query = [G.graph['canonical_map'][q] for q in query]
    G = add_supernode(G, canonical_query)
    if P is None:
        if with_restarts:
            P = compute_personalized_transition_matrix(G, alpha)
        else:
            P = compute_transition_matrix(G)
    team_ind = [G.graph['canonical_map'][c] for c in team]
    non_absorbing_nodes = [v - 1 for v in G if v not in team_ind]
    P_non_absorbing = P[non_absorbing_nodes, :][:, non_absorbing_nodes]

    step = epsilon
    s = zeros((G.number_of_nodes(), 1))
    s[-1] = 1
    X_current = s[non_absorbing_nodes]
    # Don't count the first jump, from the supernode to the queries
    score = 0
    i = 1
    while step >= epsilon :
        X_current = P_non_absorbing.T.dot(X_current)
        # Disregard the steps from the supernode to the queries
        step = X_current.sum() - X_current[-1, 0]
        score += step
        i += 1
        if max_iterations is not None and i > max_iterations:
            break
    return score
