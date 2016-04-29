## ReadMe ##

This file contains instructions on how to fill out the Execution Script for the abif_fasta function

The file "function-abif_function.R" must be downloaded and saved to the same  location as this file in order for this script to work.
Note: you do not need save these two R scripts in the same place as your ab1 files.

There are three main steps in this execution script:
 
Step 1: Set your working directory-------------------------------------------------
This tells R where to you are on the computer. Therefore, set the working directory to where you have saved this ReadMe file the execution file and the function-abif_fasta.R file. 
- Replace the square brackets with your path, in quotes; folders are separated by a forward slash (/). Example: 'C:/Users/Sandi/My Documents/'

Step 2: Sourcing required functions-------------------------------------------------
This tells R where the abif_fasta function being used is defined
- No action required. Simply run this code.

Step 3: Use the abif_fasta function------------------------------------------------
This function will read abif (ab1) files and convert it into one fasta file. It can also trim the N's from the beginning and end of the sequence.
- Fill in function's arguments according to your needs. See instructions below:

The function takes the following arguments:

folder = name of folder to where ab1 files are kept; need full path; 
         if path uses backslashes, must use double backslashes (\\) instead;
         must be in quotes.

exclude = name of files that you wish to exclude from fasta file; 
          must be in quotes, extention required, must be in a vector; 
          example: c('file1.ab1','file2.ab1')

trim = logical, default TRUE; do you want to trim N's from beginning and end?
       trim determined by looking for 4 (or other trim.level number) N's in a row.
       The end of the left trim chosen when no 4 consecutive N's detected
       The beginning of the right trim chosen 4 consecutive N's detected

trim.level = default 5; minimum number of consecutive N's to be used to trim ends;

trim.check = logical; default FALSE; print sequence before and after trim to
             manually check the trim

export.check = logical, default FALSE; export sequences before and after trimming
               for convenience when examining trim job; exports as [ouput]_check.csv

show.prog = logical, default TRUE; Reports a message on what file is currently
            being converted

ouput = default is 'V3-V6seq.FASTA'; name of output file; 
        fasta extension required; will be a FASTA file; 
        will save in location of folder provided