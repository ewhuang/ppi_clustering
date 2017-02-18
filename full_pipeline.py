### Author: Edward Huang

from multiprocessing import Pool
import os
import subprocess
import sys
import time

### Run time: About 20 minutes.

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s <species> <n_clusters>' % sys.argv[0]
        exit()
    species, n_clusters = sys.argv[1:]
    assert species in ['human'] and n_clusters.isdigit()

    command = 'build_ppi_network.py'
    print command
    subprocess.call(['python', command, species])

    # Multiple processes for simulated annealing.
    command = 'simulated_annealing.py'
    print command
    pool = Pool(processes=2) # TODO: without prosnet, only need 2 processes.
    for net_type in ['none', 'msigdb']:
        arg = ['python', command, species, net_type, n_clusters]
        pool.apply_async(subprocess.call, (arg,))
    pool.close()
    pool.join()

    # Evaluate clusters for all clusters.
    command = 'density_analysis.py'
    print command
    subprocess.call(['python', command, species, n_clusters])

    # Don't pool this for loop, since we are pooling in the script already.
    command = 'compute_label_enrichments.py'
    print command
    for label_type in ['dbgap']:
        subprocess.call(['python', command, species, n_clusters, label_type])

    # Creating cluster summarization tables.
    command = 'summarize_results.py'
    print command
    subprocess.call(['python', command, species, n_clusters])

    command = 'plot_clusters.py'
    print command
    for label_type in ['dbgap']:
        subprocess.call(['python', command, species, n_clusters, label_type])

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("Full pipeline took %s seconds..." % (time.time() - start_time))