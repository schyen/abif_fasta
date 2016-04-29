#////////////////////////////////////////////////////////////////////////////////////
# Function to read Applied Biosystem Inc. Format (Ab1, ABIF) files and put into one
# FASTA files for easy blasting on NCBI

# This script defines the function only. To execute this function,
# use the file "Execution script-abif_fasta.R"

#input:
# folder = name of folder to where ab1 files are kept; need full path; 
#          if path uses backslashes, must use double backslashes (\\) instead;
#          must be in quotes.
# exclude = name of files that you wish to exclude from fasta file; 
#           must be in quotes, extention required, must be in a vector; 
#           example: c('file1.ab1','file2.ab1')
# trim = logical, default TRUE; do you want to trim N's from beginning and end?
#        trim determined by looking for 4 (or other trim.level number) N's in a row.
#        The end of the left trim chosen when no 4 consecutive N's detected
#        The beginning of the right trim chosen 4 consecutive N's detected
# trim.level = default 4; minimum number of consecutive N's to be used to trim ends;
# trim.check = logical; default FALSE; print sequence before and after trim to
#              manually check the trim
# export.check = logical, default FALSE; export sequences before and after trimming
#               for convenience when examining trim job; exports as [ouput]_check.csv
# show.prog = logical, default TRUE; Reports a message on what file is currently
#             being converted
# ouput = default is 'V3-V6seq.FASTA'; name of output file; 
#         fasta extension required; will be a FASTA file; 
#         will save in location of folder provided

#returns a fasta file
#------------------------------------------------------------------------------------
abif_fasta <- function(folder=NULL, exclude=NULL,
                       trim=TRUE, trim.level=4, trim.check=FALSE,
                       export.check=FALSE, show.prog=TRUE, output='V3-V6seq.FASTA') {
  
  # function requirements------------------------------------------------------------
  #checking for required packages; installing those not yet installed
  if(require(seqinr)==FALSE) install.packages('seqinr')
  if(require(stringr)==FALSE) install.packages('stringr')
  
  #loading required packages
  library(seqinr)
  library(stringr)
  
  # Reading files--------------------------------------------------------------------
  #make sure folder defined
  if(is.null(folder)) stop('Name of folder not supplied. Please check the execution script preamble to make sure it is entered correctly', call.=FALSE)
  
  #making sure folder path is in correct syntax
  path <- str_replace_all(folder, '\\\\', '/')
  
  #getting abif files
  fname <- list.files(path)
  pattern <- 'ab1'
  abif_files <- str_subset(fname, pattern)
  
  #excluding specified files
  if(!is.null(exclude)) {
    
    # first check if files specified in exclude are in abif_files
    if(any(!exclude %in% abif_files)) {
      absent <- exclude[!exclude %in% abif_files]
      
      msg <- sprintf('%s could not be excluded because they were not found in your input folder. Only files in the input folder can be excluded.', 
                     str_c(absent, collapse=', '))
      stop(msg, call.=FALSE)
    }
    else{
      abif_files <- abif_files[!abif_files %in% exclude]
      msg <- sprintf('%s have been excluded from your fasta file.',
                     str_c(exclude, collapse=', '))
      message(msg)
    }
  }
  
  #reading abif files----------------------------------------------------------------
  #initializing objects
  dseq <- list()
  checkseq <- c()
  
  #going through one file at a time
  for(i in 1:length(abif_files)){
    
    fpath <- file.path(path, abif_files[i])
    abif <- read.abif(fpath)
    rawseq <- as.character(abif$Data['PBAS.1'])
    
    if(show.prog==TRUE){
      #print a message to mark progress
      msg <- sprintf('Converting file "%s"', abif_files[i])
      message(msg)
    }
    
    #Trimming------------------------------------------------------------------------
    #trimming Ns at beginning and end of sequence
    if(trim==TRUE){
      
      ind_beg <- NULL
      ind_end <- nchar(rawseq)
      
      # obtaining position of beginning trim
      for(j in trim.level:nchar(rawseq)) {
        
        # looking at each nt from the beginning and the previous 3 nt's before it
        prevnt <- c()
        for(k in rev(seq_len(trim.level-1))) {
          nt <- str_sub(rawseq, j-k, j-k)
          entry <- data.frame(prevnt=paste0('prev',k), nt=nt)
          prevnt <- rbind(prevnt, entry)
        }
        
        curr <- str_sub(rawseq, j, j)
        curr <- data.frame(prevnt='curr', nt=curr)
        prevnt <- rbind(prevnt, curr)
        
        # check if all four nt are Ns; 
        beg <- all(prevnt$nt != 'N')
        
        # if all are not Ns, then marks the beginning of seq
        if(beg == TRUE & is.null(ind_beg)) {
          ind_beg <- j-(trim.level-1)
          break
        }
      }
      
      # obtaining position of end trim
      Ns <- 0
      nt <- 0
      ind_end <- NULL
      
      # set recentN to the last N in rawseq
      lastN <- str_locate_all(rawseq, 'N')
      lastN <- as.data.frame(lastN)
      lastN <- max(lastN$end)
      
      recentN <- lastN
      
      for(m in nchar(rawseq):1) {
        # evaluate each nt, starting from the end (if is N or not N)
        curr <- str_sub(rawseq, m, m)

        # stop when number of non-Ns exceed number of Ns
        if(curr == 'N') {
          Ns <- Ns + 1
          recentN <- m
          }
        else { 
          nt <- nt + 1
          }
        
        # trim at the last N encountered
        if(nt > (Ns*2)) {
          ind_end <- recentN-1
          
          if(m < nchar(rawseq)-2) break
        }  
      
      }
      
      # trimming sequence------------------------------------------------------------
      nttrim <- str_sub(rawseq, ind_beg, ind_end)
      
      # trim quality checks and warnings---------------------------------------------
      # Checking if trim was excessive (indicates issue with sequence)
      # returns a warning
      trim_len <- nchar(nttrim)
      
      if(length(trim_len) == 0) {
        nttrim <- 'no sequence detected after trimming'
        
        msg <- sprintf('\nThe file "%s" has no remaining nucleotides after trimming.\n   Manual check of sequence is recommended!', abif_files[i])
        warning(msg, call.=FALSE)
      }
      
      else if(trim_len < 200 & length(trim_len) != 0) {
        msg <- sprintf('The file "%s" has less than 200 nucleotides after trimming.\n   Manual check of sequence is recommended!', abif_files[i])
        warning(msg, call.=FALSE)
      }
      # check if long string of Ns in sequence; return warning
      locateN <- unlist(str_locate_all(nttrim, 'N'))
      countN <- length(locateN)
      consec <- rle(diff(locateN))
      runN <- any(consec$lengths>=4 & consec$values==1)
      
      if(countN > 10 & runN == TRUE) {
        msg <- sprintf('The file "%s" has a string of Ns longer than %s after the trim. Manual check of sequence is recommended. Consider adjusting the trim.level', abif_files[i], trim.level)
        warning(msg, call.=FALSE)
      }
      
      # Reporting trim check-----------------------------------------------------------
      # print sequence and trimmed sequence for trim.check
      if(trim.check==TRUE & show.prog==FALSE) {
        msg <- sprintf('File: "%s" \n', abif_files[i])
        cat(msg)
      }
      if(trim.check==TRUE) {
        msg <- sprintf('Raw sequence:\n%s, \n\nTrimmed sequence:\n%s\n\n', 
                       rawseq, nttrim)
        cat(msg)
      }
      
      # building trim check for export-------------------------------------------------
      if(export.check==TRUE){
        entry <- data.frame(file=abif_files[i], seq_state=c('raw','trimmed'),
                            sequence=c(rawseq, nttrim))
        
        checkseq <- rbind(checkseq, entry)
      }
      
      ntseq <- nttrim
    }
    
    if(trim==FALSE){
      ntseq <- rawseq
      
      if(trim.check==TRUE) {
        message('Trim feature not in use. Set trim=TRUE to enable this feature')}
      if(export.check==TRUE) {
        message('Trim feature not in use. Set trim=TRUE to enable this feature')}
    }
    
    #putting all sequences into one object (as a list)
    dseq[[i]] <- ntseq
  }
  
  #building output file--------------------------------------------------------------
  fname_out <- file.path(path, output)
  
  if(export.check==TRUE) {
    fexport <- str_extract(output, '.*(?=\\.)')
    fexport <- paste0(fexport, '_check.csv')
    fexport <- file.path(path, fexport)
    write.csv(file=fexport, checkseq)
  }
  
  #export as one FASTA file
  write.fasta(sequences=dseq, names=abif_files, as.string=TRUE, nbchar = 1000,
              file.out=fname_out)
}
                                                                                                