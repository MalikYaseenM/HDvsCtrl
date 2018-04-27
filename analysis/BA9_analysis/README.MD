Symptomatic BA9 analysis

Mar 2018
- Find DE genes from HD BA9 and C BA9 (exclude CAP or hdpos samples) by firth
- Use top 250 genes by ascending p-adj values

Mar 29 2018

- Run fgsea on the DE genes
- Heatmap of top 250 DE genes using the unfiltered counts, cluster by samples

April 5 2018

- Find all hox genes in the hox family from the top 250 genes
- Heatmap on all the genes from the hox family only, cluster by samples
- Random forest on symptomatic BA9 top 250 DE genes
- Predict asymptomatic samples

April 12 2018

- Found 7 ribo-depleted samples from PCA (poly-A and ribo-depleted separate cleanly)
- Remove the 7 ribo-depleted samples and re-run the analysis
