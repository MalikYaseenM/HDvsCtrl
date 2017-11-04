import pandas as pd

df = pd.read_csv("clean_filtered_genes.csv")
df2 = pd.read_csv("clean_filtered_genes.csv")
total_rows = len(df.index)
start = 10499
end = 10699

# creates the sample names for snakefile to read
name_list = ""
while(total_rows > 200):
    total_rows = total_rows - 200
    name = "line" + str(start) + "to" + str(end) +" "
    start += 201
    end += 201
    name_list += name

end -= 200
name_list += "line" + str(end) + "to_end"
name_list = name_list.split()
