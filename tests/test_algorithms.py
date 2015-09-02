from __future__ import division

from networkx import is_isomorphic
from networkx.generators import karate_club_graph, krackhardt_kite_graph

from absorbing_centrality import (absorbing_centrality, add_supernode,
                                  canonical_relabel_nodes,
                                  compute_fundamental_matrix,
                                  compute_transition_matrix,
                                  greedy_team, is_canonical,
                                  keep_largest_component,
                                  update_fundamental_matrix)


def test_greedy_team():
    karate = karate_club_graph()
    solutions = {
        1: [33],
        3: [33, 0, 32],
        5: [33, 0, 32, 2, 6],
        8: [33, 0, 32, 2, 6, 25, 1, 5],
        9: [33, 0, 32, 2, 6, 25, 1, 5, 24],
        10: [33, 0, 32, 2, 6, 25, 1, 5, 24, 3]
    }
    for k in solutions:
        score_k, team_k, times_k = greedy_team(karate, k, return_times=True)
        assert sorted(team_k[-1]) == sorted(solutions[k])
