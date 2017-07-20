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
    f = csv.reader(tsvfile)
    samples_list = []
    for r in f:
        if r[0] not in samples_to_remove:
            line = str(r[0])
            split = line.split('.')
            split1 = line.split('-')
            split2 = split1[1].split('_')
            key = split2[0]
            id_list.append(key)
            samples_list.append(split[0])

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
        

