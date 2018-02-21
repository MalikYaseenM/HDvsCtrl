#!/usr/bin/env python
import pandas as pd
import sys
import os

def info_design(df, samples):
    # Takes the desired columns only, Dataset id, Subject type and subject death
    sample_info = samples.iloc[:,[14,1,4]].copy()
    sample_info.columns = ["Data_id","Subject_type","Subject_death"]

    # Makes list of hdpos ids
    hdpos = sample_info['Data_id'][sample_info['Subject_type'].str.contains("HDpos")].tolist()

    # Changes the Subject.subject_type to BA9 or CAP
    sample_info.loc[:,'Subject_type'] = sample_info["Data_id"].map(lambda x: "BA9" if "BA9" in x else "CAP")

    # Make new column list
    col = list(df)
    cols = col[:1] + hdpos

    # Drop rows if column names not in batch_id
    df.drop([col for col in df.columns if col not in cols],axis=1,inplace=True)

    # Creates new counts file with desired columns
    df.to_csv(os.path.abspath("../../samples/asymp_norm.csv"), index=False)

    ######################### Design info  #######################
    sample_i = pd.DataFrame(df.columns)
    sample_i = sample_i.drop(0)
    sample_i = sample_i.rename(columns = {0:'Data_id'})
    idesign = pd.merge(sample_i,sample_info, on='Data_id')
    # Create info design file
    idesign.to_csv(os.path.abspath("../../samples/asymp_info_design.csv"), index=False)

def main():
    if len(sys.argv) < 2:
        print("Usage: {0} <CSV File to be tab-delimited> <info-design>".format(sys.arv[0]))
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    # Read file
    df = pd.read_csv(file1, sep=',', comment='#')
    samples = pd.read_csv(file2, sep=",", comment='#')
    info_design(df, samples)

if __name__ == "__main__":
    main()
