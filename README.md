# ppi_clustering
Author: Edward Huang

## Data Pre-processing

1.
    ```bash
    python generate_directories.py
    ```

2.  Download PPI network. Go to [STRING](http://string-db.org/). Click Download, select an organism (Homo Sapiens), Update. Download the protein network data. Human data is titled 9606.protein.links.v10.txt. Move to ./data/ 

3.  Download MSigDB oncogenic gene sets. Go to [MSigDB](http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C6). Click Download entrez gene ids. Obtain c6.all.v5.2.entrez.gmt. Move to ./data/

2.  Download mappings from Entrez to ENSG IDs. Go to [biomart](http://www.ensembl.org/biomart/martview/). Choose Ensembl Genes 87, Human Genes, Attributes, uncheck Gene and Transcript ID, check Protein ID and EntrezGeneID, click Results, check Unique results only, and click Go. Downloaded file will be named mart_export.txt, rename to ensp_to_entrez.txt, move to ./data/

## Network generation

1.  
    ```bash
    python build_ppi_network.py <species>
    ```