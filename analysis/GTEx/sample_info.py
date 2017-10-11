import csv
GTEX_Sample = '../../../presymptomatic_hd_mrnaseq/samples/GTEx/GTEx_samples.txt'
absolute = '/restricted/projectnb/mlpd/DBGAP_datasets/dbGaP-14609/fastq/'
problem_samples = 'trouble_samples.txt'
out_file_name = 'GTEX_info.tsv'
stored_samples = []
samples_to_remove = []
bad_samples = []
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

output = csv.writer(open(out_file_name, 'w'), delimiter='\t',
                    quoting=csv.QUOTE_NONE, escapechar=' ')
with open(GTEX_Sample, 'r') as f:
    reader = csv.reader(f)
    for r in reader:
        if r[0] in samples_to_remove:
            pass
        else:
            line = r[0]
            path = absolute + line
            split = line.split('_')
            sample_string = split[0]+'_'+split[1]+'_'+split[2]
            if sample_string in stored_samples:
                pass
            else:
                stored_samples.append(sample_string)
                brain_region = split[1]
                sample_name = split[0].split('-')[1]
                output.writerow([sample_string, sample_name, brain_region, path])
