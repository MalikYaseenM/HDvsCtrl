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

