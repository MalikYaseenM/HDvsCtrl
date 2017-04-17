import os
import csv

fn = os.path.abspath('sample_info.csv')

ordered = open("sample_info.csv","r")
f = ordered.readlines()

for line in f:
    c = line.split()
    print(c)
ordered.close()

print(" ")

with open('sample_info_py.csv', 'r') as e:
    reader = csv.reader(e)
    for row in reader:
        print(row)
