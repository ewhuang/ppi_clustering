### Author: Edward Huang

import file_operations
import math
import sys
import matplotlib
import numpy as np
import os

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

def get_auc(point_list):
    '''
    Gets the AUC curve for a list of points. For each threshold of p-values (x-
    axis), y-axis is the number of points that has a lower p-value.
    '''
    auc_point_list = []
    # Get the sorted p-values in our list of points.
    sorted_point_list = sorted(point_list, key=lambda x: x[0])
    num_point_list = len(point_list)
    for i, (p_value, density) in enumerate(sorted_point_list):
        num_better = num_point_list - i
        auc_point_list += [(p_value, num_better)]
    return auc_point_list

def plot_clusters(plot_type):

    colors = ['blue', 'red', 'black', 'orange', 'green']
    net_type_marker_list = ('<', 'v', 'o', '^', '>')
    max_p, max_y = 0, 0

    for net_type_idx, net_type in enumerate(['none', 'msigdb']):
        point_list = []
        fname = './results/%s_results/%s_clusters_summary_%s.tsv' % (species,
            n_clusters, net_type)

        f = open(fname, 'r')
        for i, line in enumerate(f):
            line = line.strip().split('\t')
            if i == 0:
                num_gene_idx = line.index('Gene size')
                dbgap_enrichment_idx = line.index('Top DBGAP p-value')
                ratio_idx = line.index('In/(In+Out)')
                continue
            # # Skip the clusters with fewer than 30 genes. TODO:
            # if int(line[num_gene_idx]) < 30:
            #     continue

            if 'dbgap' in plot_type:
                enrichment_p = line[dbgap_enrichment_idx]
            top_enrichment_p = -math.log(float(enrichment_p), 10)
            point_list += [(top_enrichment_p, float(line[ratio_idx]))]
        f.close()

        # Plot AUC points.
        if 'auc' in plot_type:
            point_list = get_auc(point_list)

        # Update the largest y-axis value.
        if len(point_list) == 0:
            continue
        max_p = max(max_p, max(pt[0] for pt in point_list))
        max_y = max(max_y, max(pt[1] for pt in point_list))
        num_high_points = len([p for p in point_list if p[0] >= 10])
        net_type_label = '%s, %d' % (net_type, num_high_points)

        if 'auc' in plot_type:
            plt.plot(*zip(*point_list), color=colors[net_type_idx],
                label=net_type_label)
        else:
            plt.scatter(*zip(*point_list), color=colors[net_type_idx],
                label=net_type_label, marker=net_type_marker_list[net_type_idx])

    # Construct the plot details.
    plt.axvline(10)
    if 'auc' not in plot_type:
        plt.title('Cluster in-densities vs. Best enrichment p-values')
        plt.ylabel('In-density/(in-density + out-density)')
    else:
        plt.title('Number of clusters vs. Best enrichment p-values')
        plt.ylabel('Number of clusters')
    plt.xlabel('Negative log of lowest GO enrichment p-value')
    
    if 'auc' in plot_type:
        plt.legend(loc='upper right')
    else:
        plt.legend(loc='lower right')
    plt.ylim(0, max_y * 1.1)
    plt.xlim(0, max_p * 1.2)
    plt.show()

    pylab.savefig('./results/%s_results/%s_clusters_%s_plot.png' % (species,
        n_clusters, plot_type))
    plt.close()

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage: %s <species> <n_clusters> <annotation_set>' % sys.argv[0]
        exit()
    species, n_clusters, label_type = sys.argv[1:]
    assert species in ['human'] and n_clusters.isdigit()
    assert label_type in ['dbgap']

    plot_clusters(label_type)
    plot_clusters(label_type + '_auc')