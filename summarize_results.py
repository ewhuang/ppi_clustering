### Author: Edward Huang

import file_operations
import os
import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiments.

def get_cluster_densities(eval_fname):
    '''
    Get the output files of the Perl script evaluate_clustering.pl and find
    the in-density and out-density of each cluster.
    Key: Cluster Index
    Value: (in-density, out-density)
    '''
    density_dct = {}
    f = open(eval_fname, 'r')
    for i, line in enumerate(f):
        if line[:7] != 'Cluster' or 'Cluster modularity' not in line:
            continue
        line = line.split()
        clus_id, in_density, out_density = line[1], line[7], line[9]
        density_dct[clus_id] = (float(in_density), float(out_density))
    f.close()
    return density_dct

def get_enrichment_dct(enrichment_fname):
    '''
    Find the best p-value GO enrichments for each cluster.
    '''
    enrichment_dct = {}
    f = open(enrichment_fname, 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        line = line.split()
        if line[0] == 'Cluster':
            cid = line[1]
            # Skip two lines, and read in the p-values and terms.
            line = f.readline()
            term_list = line.split()
            line = f.readline()
            p_list = line.split()
            enrichment_dct[cid] = (term_list, p_list)
    f.close()
    return enrichment_dct

def write_summary_tables(net_type):
    '''
    Generates the filenames for each of the files we read from to summarize 
    the clustering performances of a particular run.
    '''
    results_folder = './results/%s_results' % species

    # Generate the format string.
    format_str = (results_folder, n_clusters, net_type)
    # Simulated annealing cluster results filename.
    if net_type == 'none':
        clus_fname = '%s/%s_clusters_%s.txt' % format_str
    else:
        clus_fname = '%s/%s_clusters_%s_clean.txt' % format_str
    clus_dct = file_operations.get_cluster_dictionary(clus_fname)
    # Perl script evaluation results filename.
    eval_fname = '%s/%s_cluster_eval_%s.txt' % format_str
    density_dct = get_cluster_densities(eval_fname)

    # DBGAP enrichment results filename.
    dbgap_fname = '%s/%s_clusters_dbgap_terms_%s.txt' % format_str
    dbgap_dct = get_enrichment_dct(dbgap_fname)

    out_fname = '%s/%s_clusters_summary_%s.tsv' % format_str
    out = open(out_fname, 'w')
    
    out.write('Cluster ID\tIn-Density\tOut-Density\tIn/(In+Out)\t'
        'Gene size\tTop DBGAP p-value\tTop DBGAP term\n')

    for cid in sorted(density_dct.keys(), key=lambda x: int(x)):
        clus = clus_dct[cid]
        num_genes = len(clus)
        in_dens, out_dens = density_dct[cid]

        best_dbgap_p = dbgap_dct[cid][1][0]
        best_dbgap_term = dbgap_dct[cid][0][0]

        out.write('%s\t%g\t%g\t%g\t%d\t%s\t%s\n' % (cid, in_dens, out_dens,
            in_dens / (in_dens + out_dens), num_genes, best_dbgap_p,
            best_dbgap_term))
    out.close()

def main():
    if len(sys.argv) != 3:
        print ('Usage:python %s <species> <n_clusters>' % sys.argv[0])
        exit()
    global species, n_clusters
    species, n_clusters = sys.argv[1:]
    assert species in ['human'] and n_clusters.isdigit()

    for net_type in ['none', 'msigdb']:
        write_summary_tables(net_type)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))