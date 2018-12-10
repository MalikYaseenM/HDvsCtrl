conda install -c conda-forge -c bioconda -f -c bubhub r python=3.5 star fastqc \
    trimmomatic multiqc snakemake pandas salmon samtools verse deeptools
conda list --export > conda_packages.txt
pip install networkx matplotlib-venn
