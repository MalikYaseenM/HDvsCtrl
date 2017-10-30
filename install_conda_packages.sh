conda install -c conda-forge -c bioconda -f -c bubhub r python=3.5 star fastqc trimmomatic=0.3.5 multiqc snakemake pandas salmon samtools verse
conda list --export > conda_packages.txt
