"""
    algorithms
    ~~~~~~~~~~

    Provides the following algorithms to pick a team of high
    absorbing centrality:
    * greedy_team
        Implements the greedy algorithm.

    :copyright: 2015 Charalampos Mavroforakis, <cmav@bu.edu> and contributors
    :license: ISC

"""


from __future__ import division

from time import clock

from networkx.algorithms.link_analysis import pagerank_scipy
from numpy import arange, zeros

from .graph_helpers import (add_supernode, canonical_relabel_nodes,
                            compute_fundamental_matrix,
                            compute_personalized_transition_matrix,
                            compute_transition_matrix, is_canonical,
                            keep_largest_component)


def greedy_team(G, k, query=None, candidates=None, fast_select=False,
                return_times=False, with_restarts=False, alpha=0.85):
    """Selects a team of nodes according to the greedy algorithm.


    Parameters
    ----------
    G : Networkx graph
        The graph from which the team will be selected.

    k : int
        The size of the team.

    query : list, optional
        If provided, the distance is measured with respect to the nodes in
        `query`.

    candidates : list, optional
        If provided, the team is picked only among the nodes in `candidates`.

    fast_select : bool, optional
        If True, the greedy algorithm will only consider candidates in a smart
        way, by examining their gain at each round (default is False).

    with_restarts : bool, optional
        If True, the greedy algorithm is based on the transition matrix w/
        restarts to the supernode (default is False).

    alpha : float, optional
        If the transition matrix has restarts, `alpha` is the probability for
        the random surfer to continue (default is 0.85).


    Returns
    -------
    scores : list
        The scores of all the greedy teams of size up to `k`.

    teams : list
        The list of greedy times of size up to `k`.

    times : list
        The time to compute each team. Returned only if `return_times` is True.
    """
    best_scores = zeros((k, 1))
    best_solutions = []
    times = []
    solution_set = []

    if not is_canonical(G):
        G = canonical_relabel_nodes(keep_largest_component(G))
    if candidates is None or not candidates:
        candidates = list(G.nodes())
    else:
        candidates = [G.graph['canonical_map'][c] for c in candidates]
    if query is None or not query:
        query = G.nodes()
    else:
        query = [G.graph['canonical_map'][q] for q in query]
    #####################################################
    # find the pagerank of the nodes before adding the supernode
    personalization = None
    if with_restarts:
        personalization = {v: 1 if v in query
                           else 0
                           for v in G}
    pagerank_centrality = pagerank_scipy(G, personalization=personalization)
    G = add_supernode(G, query)
    n = G.number_of_nodes()
    k = min(k, len(candidates))

    round_start = clock()
    if with_restarts:
        P = compute_personalized_transition_matrix(G, alpha)
    else:
        P = compute_transition_matrix(G)

    # The first node added to the solution will be the node with highest pagerank score #######
    # We start by finding the first node of the team. Instead of doing n inversions
    # we will select the node with the highest pagerank.
    current_round = 0
    candidate_values = [(c, pagerank_centrality[c]) for c in candidates]
    pagerank_candidates = [c for c, value in sorted(candidate_values,
                           key=lambda tup: tup[1], reverse=True)[:5]]
    round_min = -1
    round_best = pagerank_candidates[0]
    for pagerank_candidate in pagerank_candidates:
        non_absorbing_nodes = [i for i in arange(n)
                               if i + 1 not in [pagerank_candidate]]
        P_abs = P[non_absorbing_nodes, :][:, non_absorbing_nodes]
        F = compute_fundamental_matrix(P_abs, fast=True)
        row_sums = F.sum(axis=1)
        score = row_sums[-1].sum() - F[-1, -1]
        if score < round_min or round_min == -1:
            round_min = score
            round_best = pagerank_candidate
    best_candidate = round_best
    solution_set.append(G.graph['label_map'][best_candidate])
    times.append(clock() - round_start)

    non_absorbing_nodes = [i for i in arange(n)
                           if i + 1 != best_candidate]
    P_abs = P[non_absorbing_nodes, :][:, non_absorbing_nodes]
    F = compute_fundamental_matrix(P_abs, fast=True)

    best_scores[current_round] = round_min
    best_solutions.append(solution_set)
    current_round += 1
    # We will take advantage of the submodularity of the problem
    # At each round, we remember ,for ech node, the gain that adding this node has
    # At the next round, we compute a running round_best_gain (e.g. from the first node in the loop)
    # and for every next node, we first check if that node's gain in the previous round was smaller
    # than the current round_best_gain. The submodularity property suggests that if this is the case,
    # then in the previous values is an upper bound for this round's gain for that node, so we skip checking it
    round_gain = {}
    nodes_to_check = candidates
    absorbing_nodes = [best_candidate]
    while len(solution_set) < k:
        round_min = -1
        round_best_member = -1
        round_start = clock()
        for c in nodes_to_check:
            if G.graph['label_map'][c] in solution_set:
                continue
            if fast_select and c in round_gain:
                gain_upper_bound = round_gain[c]
                if gain_upper_bound < best_scores[current_round - 1] - round_min:
                    # the node's gain cannot exceed its upper bound
                    continue
            non_absorbing_nodes = [i for i in arange(n)
                                   if i + 1 not in absorbing_nodes + [c]]
            P_abs = P[non_absorbing_nodes, :][:, non_absorbing_nodes]
            F_new = compute_fundamental_matrix(P_abs)
            # Now that we inverted P, we need to compute the absorbing time.
            row_sums = F_new.sum(axis=1)
            total_steps = row_sums[-1].sum() - F_new[-1, -1]
            round_gain[c] = best_scores[current_round - 1] - total_steps
            if total_steps < round_min or round_min == -1:
                round_best_member = c
                round_min = total_steps
        absorbing_nodes.append(round_best_member)
        solution_set = solution_set + [G.graph['label_map'][round_best_member]]
        best_scores[current_round] = round_min
        best_solutions.append(solution_set)
        if fast_select:
            # order the candidates by decreasing gain,
            # so that we increase the chances to find the next node faster
            nodes_to_check = sorted(round_gain,
                                    key=round_gain.get,
                                    reverse=True)
        times.append(clock() - round_start)
        current_round = current_round + 1
    if return_times:
        return (best_scores, best_solutions, times)
    else:
        return (best_scores, best_solutions)
