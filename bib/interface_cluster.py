import ctl

import community
import networkx as nx


'''
Module providing functions for clustering of interface matrices.
'''

def cluster(corr_coefficients, mode=0, verbous=False):
    '''
    Clustering according to correlation coefficients.
    '''
    g = nx.MultiGraph()
    nodes = []
    edges = []

    # proportionality constant between edge weight and correlation coefficient
    alpha = 1

    for c in corr_coefficients:
        rank_id0 = int(c[0][3]. \
                       replace('unrelaxed_rank', '').replace('rank', ''). \
                       split('_')[0])
        rank_id1 = int(c[1][3]. \
                       replace('unrelaxed_rank', '').replace('rank', ''). \
                       split('_')[0])

        # define node ids
        n0 = c[0][0].split('/')[-1]+'_'+str(c[0][2])+'_'+str(rank_id0)
        n1 = c[1][0].split('/')[-1]+'_'+str(c[1][2])+'_'+str(rank_id1)

        if n0 not in nodes:
            nodes.append(n0)

        if n1 not in nodes:
            nodes.append(n1)

        minarea = 0

        if corr_coefficients[c][0] > 0 and \
           corr_coefficients[c][1] >= minarea:

            # modes:
            # 0: set edge weights to correlation coeff
            # 1: set all edge weights to 1
            if mode == 0:
                edge = (min(n0, n1), max(n0, n1), alpha*corr_coefficients[c][0])
            elif mode == 1:
                edge = (min(n0, n1), max(n0, n1), 1)
            else:
                ctl.error('cluster: mode not set correctly')

            if edge not in edges:
                edges.append(edge)

    for n in nodes:
        g.add_node(n)

    for e in edges:
        g.add_edge(e[0], e[1], weight=e[2])

    if verbous:
        ctl.d('edge weights:')
        for e in g.edges(data=True, keys=True):
            ctl.d(e)

    partitions = community.best_partition(g, weight='weight', resolution=1)
    partitions_n = 0

    for p in partitions:
        if partitions[p] > partitions_n:
            partitions_n = partitions[p]

    partitions_n += 1

    return partitions, g, partitions_n


def sort_clusters(partitions, partitions_n):
    ''' Sort clusters by the highest score that occurs for a node. '''

    partitions_ = []
    partitions_maxscore = [[-1, i] for i in range(0, partitions_n)]
    partitions_replacementtable = [-1 for i in range(0, partitions_n)]

    for p in partitions:
        score = float(p.split('_')[-2])
        if score < 0 or score > 1:
            ctl.error('sort_cluster: score not in range')

        partitions_maxscore[partitions[p]][0] = \
                            max(partitions_maxscore[partitions[p]][0], score)
            
    partitions_maxscore_sorted = sorted(partitions_maxscore, \
                                        key=lambda x: x[0], reverse=True)

    for i,p in enumerate(partitions_maxscore_sorted):
        partitions_replacementtable[p[1]] = i

    partitions_sorted = {}

    for p in partitions:
        partitions_sorted[p] = partitions_replacementtable[partitions[p]]

    return partitions_sorted


def plot_clusters(g, fn, plt):
    '''
    Plot clusters.
    '''
    plt.figure(figsize=(12, 8))
    layout = nx.spring_layout(g, k=2.0)
    nx.draw(g, layout, \
            font_size=8, edge_color='gray', \
            with_labels=True)

    plt.savefig(fn, dpi=600)

    return
