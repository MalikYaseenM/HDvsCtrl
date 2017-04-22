import pandas as pd

## List of bad genes so far:
# Gene name, Line from all_salmon_norm_filter (in excel), line # pandas
# For excel, first row is header, genes start from second row
#1  ENSG00000280594.1,  line 237,     ,line 235
#2  ENSG00000197405.7,  line 3279     ,line 3277
#3  ENSG00000149257.13, line 3302     ,line 3300
#4  ENSG00000113263.12, line 5571     ,line 5569
#5  ENSG00000162772.16, line 9679     ,line 9677
#6  ENSG00000174576.8,  line 10001    ,line 9999
#7  ENSG00000173110.7,  line 10498    ,line 10496
#8  ENSG00000118193.11, line 11866    ,line 11864
#9  ENSG00000163874.9,  line 12358    ,line 12356
#10 ENSG00000132002.7,  line 15059    ,line 15057
#11 ENSG00000166592.11, line 15463    ,line 15461
#12 ENSG00000275302.1,  line 15488    ,line 15486
#13 ENSG00000087074.7,  line 21133    ,line 21131
#14 ENSG00000253190.3,  line 21905    ,line 21093

df = pd.read_csv("all_salmon_norm_filter.csv")
# drops all bad gene rows only
# pandas index starts from 0
df = df.drop([235, 3277, 3300, 5569, 9677, 9999, 10496, 11864, 12356, 15057, 15461, 15486, 21131, 21903])

# Firth testing
##drops row 0 to X returns rows starting from X
##drop bad gene #1
# df = df.drop(df.index[0:236])
##drop bad gene #2
# df = df.drop(df.index[0:3278])
##drop bad gene #3
# df = df.drop(df.index[0:3301])
##drop bad gene #4
# df = df.drop(df.index[0:5570])
##drop bad gene #5
# df = df.drop(df.index[0:9678])
##drop bad gene #6
# df = df.drop(df.index[0:10000])
##drop bad gene #7
# df = df.drop(df.index[0:10497])
##drop bad gene #8
# df = df.drop(df.index[0:11865])
##drop bad gene #9
# df = df.drop(df.index[0:12357])
##drop bad gene #10
# df = df.drop(df.index[0:15058])
##drop bad gene #11
# df = df.drop(df.index[0:15462])
##drop bad gene #12
# df = df.drop(df.index[0:15487])
##drop bad gene #13
# df = df.drop(df.index[0:21132])
##drop bad gene #14
#df = df.drop(df.index[0:21904])

# write to csv file
df.to_csv("clean_filtered_genes.csv", index=False)

