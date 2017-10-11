This is where the Poly A vs Ribo-depleted analysis will took place.
# Step 1
protocol.py - Takes the file HD_mRNASeq_sample_info.csv and gets the columns, Dataset.dataset_id', and 'Dataset.protocol'.

Dataset.dataset_id contains all the sample id names.
Dataset.protocol describes what library prep method was used on that sample. The protocols were ribo-depletion and poly-A selection. 

These values are then written to a csv called sample_info.csv

# Step 2 
workflow.py - takes the newly created sample_info.csv and seperates the samples out into seperate csv based on whether they are poly-A or ribo-depleted/

# Step 3
The snakefile is what is used to actually run all of the python scrips.
protocol.py is run by the snakefile to create the sample_info.csv using the rule Protocol_Appender

The PolyRibocsv rule runs the workflow.py file which creates the two seperate poly-A and ribo-depleted csv's. 

After running the snakefile you should have two csv files in the output. polyseq.csv, and ribodepseq.csv. polyseq.csv will contain 1 column called 'Poly_A' with rows equal to the ids that were poly-A selected.

ribodepseq.csv will contain 1 column labled Ribo_D and will have the names of all of the samples that were ribo-depleted.



