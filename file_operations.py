### Author: Edward Huang

### Contains functions that read globally used files.

def read_entrez_to_ensp_dct():
    '''
    Returns the mappings from entrez to ENSP as a dictionary.
    Key: EntrezGene ID -> str
    Value: set of ensembl protein IDs (ENSP) -> set(str)
    '''
    entrez_to_ensp_dct = {}
    f = open('./data/ensp_to_entrez.txt', 'r')
    for i, line in enumerate(f):
        if i == 0: # Skip header line.
            continue
        line = line.split()
        if len(line) != 2:
            continue
        ensp_id, entrez_id = line
        assert entrez_id.isdigit()
        if entrez_id not in entrez_to_ensp_dct:
            entrez_to_ensp_dct[entrez_id] = set([])
        entrez_to_ensp_dct[entrez_id].add(ensp_id)
    f.close()
    return entrez_to_ensp_dct

def read_oncogenic_signatures(ppi_protein_set):
    '''
    Returns the oncogenic signatures provided by MSigDB. Returns a dictionary.
    Key: Cellular pathways that are often dis-regulated in cancer -> str
    Value: sets of genes -> set(str)
    '''
    entrez_to_ensp_dct = read_entrez_to_ensp_dct()
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