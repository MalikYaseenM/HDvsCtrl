This directory contains the pipeline for the processing of data Downloaded from 
https://www.gtexportal.org/home/

# General Information
This directory contains files from the GTEX database.
There are 496 source files. They are paired end so there are 248 Read1's
and 248 Read2's 
After looking at the multiqc report 112 samples were filtered out
This left 384 files (192 pairs)
A tsv file called GTEX_info.tsv contains information which can be parsed to work with the samples.



Diego Crespo originally worked on these files
Email: crespodi@bu.edu

# Donor Inclusion/Exclusion Criteria are as follows:
1. 21 ≤ Age (years) ≤ 70
2. 18.5 < BMI < 35
3. Time between death and tissue collection less than 24 hours
4. No whole blood transfusion within 48 hours prior to death
5. No metastatic cancer
6. No chemotherapy or radiation therapy within the 2 years prior to death
7. Generally unselected for presence or absence of diseases or disorders, except for potentially communicable diseases that disqualify someone to donate organs or tissues would also be disqualifying for GTEx.


#Info
The First step was to run the snakefile and process all of the GTEX samples and create a multiqc report. From the Multiqc report
trouble_samples.txt was created. It contains a list of samples which after processing all of the data from GTEX had various issues.

GTEX-T5JC_BA9_mRNASeq_R1 and GTEX-T5JC_BA9_mRNASeq_R2

GTEX-T6MN_BA9_mRNASeq_R1 and GTEX-T6MN_BA9_mRNASeq_R2
GTEX-XLM4_BA9_mRNASeq_R1 and GTEX-XLM4_BA9_mRNASeq_R2

were disqualified due to having lengths of 250bp in a set that otherwise had lengths of 76bp.


The rest of the samples were disqualified due to have qualitatively different properties from the majority of the samples in the Per Base Sequence Content. In total 53 paired end samples were disqualified for further analysis after processing.

HD_Salmon.txt contains a list of the control and HD samples that were used in a previous analysis that have not been used yet with GTEX but are there just in case.

#Snakefile

The Snakefile contains the pipeline that was used to process all of the GTEx data.

The pipeline consists of Fastqc, Trimmomatic, Salmon, STAR, and Multiqc. By the end  of this analyis, the samples that we were going to use in the differential gene expression were established.







