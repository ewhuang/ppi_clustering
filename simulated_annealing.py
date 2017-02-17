### Author: Edward Huang

import file_operations
import os
import subprocess
import sys
import time

def main():
    if len(sys.argv) != 4:
        print ('Usage:python %s <species> <net_type> <n_clusters>' % 
            sys.argv[0])
        exit()
    species, net_type, n_clusters = sys.argv[1:]
    assert species in ['human'] and net_type in ['none', 'msigdb']
    assert n_clusters.isdigit()

    binary = './WlogV/makedir/WlogVImplement'
    temp = 10 # TODO.
    command = ('%s %s 1 0 "./data/%s_data/orth.txt" 1 "./data/%s_data/'
        'network_%s.txt" -t %s 2>log > "./results/%s_results/%s_clusters_'
        '%s.txt"' % (binary, n_clusters, species, species, net_type, temp,
            species, n_clusters, net_type))
    print command
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))