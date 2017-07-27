import csv
#samples = f.open('samples.txt', 'r')
# ['GTEX', 'ZF28_BA9_mRNASeq_R1.fastq.gz']
# string1 = "GTEX-XLM4_HYP_mRNASeq_R1.fastq.gz"
# string2 = "GTEX-XLM4_NUA_mRNASeq_R1.fastq.gz"
#split = ['GTEX-11GSP_BA9_mRNASeq_R1', 'fastq', 'gz'] with '.'
# split = string1.split('-')
# strip = split[1].split('_')
#files = '../../../presymptomatic_hd_mrnaseq/samples/GTEx/four_samples.txt'
#problem_samples = 'bad_samples.txt'
path = '/usr3/graduate/crespodi/Huntington/presymptomatic_hd_mrnaseq/samples/'
problem_samples = 'trouble_samples.txt'
bad_samples = []
samples_to_remove = []
with open(problem_samples, 'r') as textfile:
    f = csv.reader(textfile)
    for r in f:
        line = str(r[0])
        sample = line[0:-3]
        bad_samples.append(sample)
processed_samples = set(bad_samples)
for samples in processed_samples:
    individual1 = samples+'_R1.fastq.gz'
    individual2 = samples+'_R2.fastq.gz'
    samples_to_remove.append(individual1)
    samples_to_remove.append(individual2)
    
files = '../../../presymptomatic_hd_mrnaseq/samples/GTEx/GTEx_samples.txt'
f = csv.reader(files)
id_list = []
with open(files, 'r') as tsvfile:
    samples = []
    f = csv.reader(tsvfile)
    for r in f:
        if r[0] not in samples_to_remove:
            line = str(r[0])
            split = line.split('.')
            individual = split[0][0:-3]
            if individual in samples:
                pass
            else:
                samples.append(individual)
                
# hd_samples = 'HD_Salmon.txt'
# total_samples = []
# with open(hd_samples, 'r') as csvfile:
#     f = csv.reader(csvfile)
#     for r in f:
#         line = path +str(r[0])+'/quant.genes.sf'
#         print(line)
        
# quant_sf = '../../../presymptomatic_hd_mrnaseq/samples/GTEx/genes_with_quantsf.txt'
# with open(quant_sf, 'r') as csvfile:
#     quant_samples = []
#     f = csv.reader(csvfile)
#     for r in f:
#         line = str(r[0])
#         split = line.split('__')
#         quant_samples.append(split[0])
# quant_set = set(quant_samples)
# for sample in samples:
#     if sample not in quant_set:

#coord_list = ['GTEX-N7MT_CAU_mRNASeq_R1', 'GTEX-WHSE_CAU_mRNASeq_R1', 'GTEX-XLM4_BA9_mRNASeq_R2', 'GTEX-XMD1_BA9_mRNASeq_R1']
# trouble_list = []
# for i in range(len(samples_list)):
#     if i in range(298, 399):
#         trouble_list.append(samples_list[i])
#     elif i in range(435,437):
#         trouble_list.append(samples_list[i])

# trouble_file = open('trouble_samples.txt', 'w')
# for i in range(len(trouble_list)):
#     trouble_file.write("%s\n" % trouble_list[i])
        

