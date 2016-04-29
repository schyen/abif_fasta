#////////////////////////////////////////////////////////////////////////////////////
# Use this script to read Applied Biosystem Inc. Format (Ab1, ABIF) files and 
# put into one FASTA files for easy blasting on NCBI.

# To execute this code, fill out the required information according to the 
# instructions in the ReadMe file.
# When ready, highlight the whole script and click "Run"

# Step 1: Set your working directory-------------------------------------------------

path <- '[your/path/here]'
setwd(path)

# Step 2: Sourcing required functions------------------------------------------------

func_path <- file.path(path, 'function-abif_fasta.R')
source('function-abif_fasta.R')

# Step 3: Use the abif_fasta function------------------------------------------------

abif_fasta(folder=NULL, exclude=NULL, trim=TRUE, trim.level=4, 
           trim.check=FALSE, export.check=FALSE,
           show.prog=TRUE, output='V3-V6seq.FASTA')

