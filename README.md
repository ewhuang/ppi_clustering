# ppi_clustering
Author: Edward Huang

## Data Pre-processing

1.  Generate directories for further pre-processing.

    ```bash
    python generate_directories.py
    ```

2.  Download PPI network. Go to [STRING](http://string-db.org/). Click Download, select an organism (Homo Sapiens), Update. Download the protein network data. Human data is titled 9606.protein.links.v10.txt. Move to ./data/ 

3.  Download MSigDB oncogenic gene sets. Go to [MSigDB](http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C6). Click Download entrez gene ids. Obtain c6.all.v5.2.entrez.gmt. Move to ./data/

2.  Download mappings from Entrez to ENSG IDs. Go to [biomart](http://www.ensembl.org/biomart/martview/). Choose Ensembl Genes 87, Human Genes, Attributes, uncheck Gene and Transcript ID, check Protein ID and EntrezGeneID, click Results, check Unique results only, and click Go. Downloaded file will be named mart_export.txt, rename to ensp_to_entrez.txt, move to ./data/

## Network generation

1.  Create the network to be clustered on by simulated annealing.

    ```bash
    python build_ppi_network.py <species>
    ```

## Clustering

1.  Compile simulated annealing with WlogV objective function code.

    ```bash
    $ cd makedir
    $ rm *
    $ cmake ..
    $ make
    ```

2.  Run simulated annealing script.

    ```bash
    python simulated_annealing.py <species> <network_type> <n_clusters>
    ```

3.  Evaluate clusters.

    ```bash
    python density_analysis.py <species> <n_clusters>
    ```