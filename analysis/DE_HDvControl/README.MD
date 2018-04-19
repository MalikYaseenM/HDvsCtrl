Nov 4 2017

HD versus Control Analysis version 2

Problem: Batch effect of poly-A and ribo depleted RNA processing.
- Plot all normalized counts in PCA and found that it separates evenly into Poly-A or ribo-depleted.
- Found an outlier poly-A sample in ribo-depleted cluster: 'H_0014_BA9_mRNASeq'

Solution: Remove 'H_0014_BA9_mRNASeq' from the analysis.

Steps for analysis (See Snakefile and python files for more info):
Snakefile pipeline only covers steps 1-3.
1. Remove H_0014 sample
2. Filter the rest of the samples
2. Normalize
3. Firth
4. Analysis of DE genes from firth 
5. DAVID

Other analysis:
- Separating poly-A and ribo depleted and analyzing them separately

Steps for poly-A only analysis (poly-A snakefile):
Snakefile pipeline only covers steps 1-4.
1. Remove H_0014 sample and ribo-depleted processed samples
2. Filter 
3. Normalize
4. Firth
5. Analysis of DE genes from firth
6. DAVID

Steps for ribo-depleted only analysis (riboD snakefile):
Snakefile pipeline only covers steps 1-4.
1. Remove all poly-A processed samples
2. Filter
3. Normalize
4. Firth
5. Analysis of DE genes from firth
6. DAVID

Feb 2018
- new pol-A samples -> replaces the ribo-depleted 
- re-run analysis

Apr 2018
- Still found 7 ribo-depleted samples
- Remove the 7 ribo-depleted samples
-  Re-run analysis

