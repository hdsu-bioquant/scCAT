# scCAT
Deconvolution of single-cell multi-omics layers reveals regulatory heterogeneity 

Code for inferring regulatory relationship


# scCATseq Regulatory Relationships
Analysis includes following steps:
1. Preprocessing
2. Non-Negative Matrix Factorization
3. Computation of active TFs using SCENIC
4. Computation of correlation between promoter peaks and distal peaks signal
5. Computation of correlation between gene expression and distal peaks signal
6. Cpmputation of cell specific Regulatory Relationships

# HOW TO USE:
1. Install SNAKEMAKE, highly recommended to follow this [TUTORIAL](http://snakemake.bitbucket.org/snakemake-tutorial.html) using Miniconda3
2. Adjust configs/config.yaml file.

    2.1. Insert your space delimated assay identifiers

    2.2. Change Path to your data and Results folder. 

    REQUIRED DATA STRUCTURE:
    
    ```bash
    INPUT:
    data/
    ├── CellLines/
    │   ├── annotation/
    │   │   └── CellLines_metadata.Rds
    │   ├── atac/
    │   │   ├── CellLines_ATACseqCounts.Rds.Rds
    │   │   └── CellLines_ATACseqGRanges.Rds
    │   └── rna/
    │       └── CellLines_RNAseqCounts.Rds
    ├── HumanEmbryo/

    OUTPUT:
    resutls/
    ├── CellLines/
    │   ├── atac/
    │   │   └──
    │   ├── rna/
    │   │   ├──
    │   │   └──
    │   └── rna_atac/
    │       └──
    ├── HumanE
    ```

    2.3. OPTIONAL: Select which rules of the graph you want to run, by changing Trua and False flags in the "Pipeline Steps to Run" section

    2.4. OPTIONAL: Change parameters used to run NMF and to infer regulatory interactions


3. Make a Snakemake try run: 
```bash
# If your are using Miniconda3
source activate <YOUR SNAKEMAKE ENV NAME>
# Perform snakemake dry run
snakemake -np 
```

4. The pipeline could be run in a cluster, make sure to modify the cluster/qsub.yaml file or create a new one for your system specification

5. The pipeline uses conda environments to install almost all dependencies. In the case of the rules using GPUs, your system should have installed the required cuda drivers.

6. Check if all the rules you want, are included in the snakemake graph:
```bash
# If your are using Miniconda3
source activate <YOUR SNAKEMAKE ENV NAME>
# Print rule graph
snakemake -F  --rulegraph|dot -Tpdf > rulegraph.pdf
```

7. If the dry run is working and everything you need will be created, run the pipeline on the cluster: 
```bash
# If your are using Miniconda3
source activate <YOUR SNAKEMAKE ENV NAME>
# Perform snakemake run (using qsub)
snakemake -k -w 50 --jobs 50 --use-conda  --cluster-config cluster/qsub.yaml --cluster "qsub -l 'walltime={cluster.walltime}, nodes=1:ppn={cluster.cores}, mem={cluster.memory}'"

```






