Feb 21 2018

New Asymptomatic analysis pipeline and scripts
- samples are now all pol-A processed

1. Make info design file
2. Get the 5 asymptomatic individual normalized counts from the all_norm.csv
3. Run firth
4. Since the sample size is so small, it's unlikely that there will be genes with significant p values
5. Use jupyter notebook to analyze the dataset from step 2

March 22 2018
- Do asymptomatic HD BA9 versus control BA9 matching by age
- Run firth, use counts beta on fgsea R
- Get the samples from all_norm.csv, filter by outliers, perform t-stat
- Compare 'DE' results from filtering method with DE GTEx Control BA9 CAP
