import pandas as pd

df = pd.read_csv("clean_filtered_genes.csv")
df2 = pd.read_csv("clean_filtered_genes.csv")
total_rows = len(df.index)
start = 0
end = 200

# 200 rows/file
while(total_rows > 200):
    df1 = df.drop(df.index[201:])
    df = df.drop(df.index[0:201])
    total_rows = total_rows - 200
    name = "line" + str(start) + "to" + str(end) + ".csv"
    start += 201
    end += 201
    df1.to_csv(name, index=False)

end -= 200
last_line = "line" + str(end) + "to_end.csv"
df.to_csv(last_line, index=False)
