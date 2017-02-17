### Author: Edward Huang

from file_operations import read_oncogenic_signatures
import os
import sys
import time

### This script takes the PPI network and creates the network files ready for
### simulated annealing. 

def generate_directories():
    global data_folder, results_folder
    data_folder = './data/%s_data' % species
    results_folder = './results/%s_results' % species
    for folder in [data_folder, results_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

def read_protein_links():
    '''
    Returns the PPI network from the STRING database as a set of tuples. Also
    returns the set of proteins in the network.
    '''
    protein_edge_set, protein_set = set([]), set([])
    # Select filename.
    if species == 'human':
        fname = '9606.protein.links.v10'

    f = open('./data/%s.txt' % fname, 'r')
    for line in f:
        if 'ENSP' not in line: # Skip header line.
            continue
        edge = line.split()[:2]

        # Get rid of leading '9606.' substrings.
        for i, protein in enumerate(edge):
            edge[i] = protein[protein.index('.')+1:]
            protein_set.add(edge[i])

        edge = tuple(edge) # Convert to tuple for hashability.
        # Skip duplicate reverse edges.
        if edge[::-1] in protein_edge_set:
            continue
        # Update edge dictionary.
        protein_edge_set.add(edge)
    f.close()
    return protein_edge_set, protein_set

def write_orth_file(gene_example):
    '''
    Writes out the orthology file.
    '''
    # Write out the orth file.
    orth_out = open('%s/orth.txt' % data_folder, 'w')
    orth_out.write('0\t0\t%s\t%s' % (gene_example, gene_example))
    orth_out.close()

def write_network(protein_edge_set, protein_set, net_type):
    '''
    Given a protein edge dictionary, write them out to prepare for clustering.
    If the network type denotes "with msigdb", then we also add in the
    oncogenic signature labels as edges. Also write out to the "real" evaluation
    network.
    '''
    num_nodes = len(protein_set)
    if net_type == 'msigdb':
        oncogenic_signature_dct = read_oncogenic_signatures(protein_set)
        num_nodes += len(oncogenic_signature_dct)

    net_out = open('%s/network_%s.txt' % (data_folder, net_type), 'w')
    net_out.write('0\n%d\n' % num_nodes)

    # "Real" evaluation network.
    real_out = open('%s/eval_network_%s.txt' % (data_folder, net_type), 'w')
    real_out.write('Real network\n')

    for (p_1, p_2) in protein_edge_set:
        # Write out the edge twice.
        net_out.write('%s\t%s\t1\n%s\t%s\t1\n' % (p_1, p_2, p_2, p_1))
        real_out.write('0\t%s\t%s\t1\n0\t%s\t%s\t1\n' % (p_1, p_2, p_2, p_1))

    if net_type == 'msigdb':
        for gsn in oncogenic_signature_dct: # gsn stands for gene set name.
            protein_set = oncogenic_signature_dct[gsn]
            for p in protein_set: # p stands for protein
                # write each lien out twice.
                net_out.write('%s\t%s\t1\n%s\t%s\t1\n' % (p, gsn, gsn, p))
                real_out.write('0\t%s\t%s\t1\n' % (p, gsn))
                real_out.write('0\t%s\t%s\t1\n' % (gsn, p))
    net_out.close()
    real_out.close()

    # Write out an aribtrary protein to the orthology file.
    write_orth_file(p_1)

def main():
    if len(sys.argv) != 2:
        print 'Usage: python %s <species>' % sys.argv[0]
        exit()
    global species
    species = sys.argv[1]
    assert species in ['human']

    generate_directories()

    protein_edge_set, protein_set = read_protein_links()

    write_network(protein_edge_set, protein_set, 'none')
    write_network(protein_edge_set, protein_set, 'msigdb')

if __name__ == '__main__':
    start_time = time.time()
    main()
    print "---%f seconds---" % (time.time() - start_time)