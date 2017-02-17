### Author: Edward Huang

from multiprocessing import Pool
import os
import sys
import time

### This script calls the perl script to evaluate clusters for a given run.
### Outputs the file to a results folder.
### Run time: 4.5 minutes

def remove_non_genes():
    '''
    Rewrite the clustering output, but only preserving nodes that are proteins.
    '''
    f = open('%s/%s_clusters_msigdb.txt' % (results_folder, n_clusters), 'r')
    out = open('%s/%s_clusters_msigdb_clean.txt' % (results_folder, n_clusters),
        'w')
    for i, line in enumerate(f):
        if i == 0 or 'ENSP' in line:
            out.write(line)
    f.close()
    out.close()

def call_perl_script(net_type):
    suffix = net_type[:]
    if net_type == 'msigdb':
        suffix = 'msigdb_clean'
    clus_fname = '%s/%s_clusters_%s.txt' % (results_folder, n_clusters, suffix)

    command = ('perl ./evaluate_clustering.pl "%s" "./data/%s_data/real_'
        'network_%s.txt" > "%s/%s_cluster_eval_%s.txt"' % (clus_fname, species,
            net_type, results_folder, n_clusters, net_type))
    print command
    os.system(command)

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s <species> <n_clusters>' % (sys.argv[0])
        exit()
    global species, n_clusters
    species, n_clusters = sys.argv[1:]
    assert species in ['human'] and n_clusters.isdigit()

    global results_folder
    results_folder = './results/%s_results' % species

    remove_non_genes()

    pool = Pool(processes=2)
    # Simulatenously run both scripts.
    pool.apply_async(call_perl_script, ('none',))
    pool.apply_async(call_perl_script, ('msigdb',))

    pool.close()
    pool.join()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))