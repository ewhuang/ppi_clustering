### Author: Edward Huang

### Contains functions that read globally used files.

def read_external_to_ensp_dct(external_id_type):
    '''
    Returns the mappings from external to ENSP as a dictionary.
    Key: externalGene ID -> str
    Value: set of ensembl protein IDs (ENSP) -> set(str)
    '''
    external_to_ensp_dct = {}
    f = open('./data/ensp_to_%s.txt' % external_id_type, 'r')
    for i, line in enumerate(f):
        if i == 0: # Skip header line.
            continue
        line = line.split()
        if len(line) != 2:
            continue
        ensp_id, external_id = line
        if external_id not in external_to_ensp_dct:
            external_to_ensp_dct[external_id] = set([])
        external_to_ensp_dct[external_id].add(ensp_id)
    f.close()
    return external_to_ensp_dct

def read_oncogenic_signatures(ppi_protein_set):
    '''
    Returns the oncogenic signatures provided by MSigDB. Returns a dictionary.
    Key: Cellular pathways that are often dis-regulated in cancer -> str
    Value: sets of genes -> set(str)
    '''
    entrez_to_ensp_dct = read_external_to_ensp_dct('entrez')
    oncogenic_signature_dct = {}
    f = open('./data/c6.all.v5.2.entrez.gmt', 'r')
    for line in f:
        line = line.split()
        gene_set_name, url, gene_set = line[0], line[1], line[2:]
        assert 'http://' in url and gene_set_name not in oncogenic_signature_dct
        # Convert from Entrez to ENSP.
        mapped_gene_set = set([])
        for entrez_id in gene_set:
            if entrez_id in entrez_to_ensp_dct:
                ensp_id_set = entrez_to_ensp_dct[entrez_id]
                mapped_gene_set = mapped_gene_set.union(ensp_id_set)
        # Only keep protein mappings that occur in the PPI network.
        mapped_gene_set = mapped_gene_set.intersection(ppi_protein_set)
        if len(mapped_gene_set) == 0:
            continue
        oncogenic_signature_dct[gene_set_name] = mapped_gene_set
    f.close()
    return oncogenic_signature_dct

def get_cluster_dictionary(filename):
    '''
    Returns a dictionary of clusters.
    Key: cluster ID -> str
    Value: lists of genes in the cluster-> list(str)
    '''
    cluster_dct = {}
    f = open(filename, 'r')
    # Read in the cluster file to create the cluster dictionary.
    for i, line in enumerate(f):
        if i == 0:
            continue
        newline = line.strip().split('\t')
        cluster = newline[2][len('Cluster '):]
        # Skip garbage clusters.
        if cluster == '0':
            continue
        gene = newline[1][len('Gene '):]
        assert 'ENSP' in gene
        if cluster not in cluster_dct:
            cluster_dct[cluster] = []
        cluster_dct[cluster] += [gene]
    f.close()
    return cluster_dct