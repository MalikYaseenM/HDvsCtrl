import pandas as pd
info = open('sample_info.csv','r')
fullread = info.readlines()

Poly = []
Ribo = []

for lines in fullread:
    if 'poly-A' in lines:
        newline= lines.strip('\n')
        split = newline.split(',')
        Poly.append(split[1])
    else:
        newline= lines.strip('\n')
        split = newline.split(',')
        Ribo.append(split[1])

polyseq=pd.DataFrame(Poly, columns=['Poly_A'])
polyseq.to_csv('polyseq.csv', sep=',', index=False)
ribodepseq=pd.DataFrame(Ribo, columns=['Ribo_D'])
ribodepseq.to_csv('ribodepseq.csv', sep=',', index=False)
