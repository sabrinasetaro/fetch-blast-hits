#R script written by Sabrina Setaro
#File to parse a blast xml file, fetch sequences from hits and save them to a fasta file

library(xml2)
library(seqinr)
library(plyr)

#############################################
# getting the data

#get all data
xml_files <- c("file1.xml", "file2.xml")

#test which ones returned no results
skipNoHits <- function (file) {
  url <- xml_files[file]
  xmlfile <- read_xml(url)
  noHits <- xml_text(xml_find_all(xmlfile, xpath = "//Iteration_message"))
  if(length(noHits)==0) {
    return(url)
  }
}

xml_with_hits <- unlist(sapply(1:length(xml_files), skipNoHits))

getDatafromXML <- function(file) {

  #get xml file
  url <- xml_with_hits[file]
  #url <- "AT1G05055.xml"
  xmlfile <- read_xml(url)
  
  #############################################
  # Functions and Helpers
  
  #function to get all sample names to find out duplicates
  splitall <- function(x) {
    split.ls <- strsplit(ids[[x]], "_")
    return(split.ls[[1]][[2]])
  }
  
  #function to get info about how many hsps are combined in hit
  getHsps <- function(x) {
    xml_length(xml_find_all(xmlfile, xpath = "//Hit_hsps")[[x]])
  }
  
  #get number of hits
  hits <- length(xml_find_all(xmlfile, xpath = "//Hit"))
  
  #get hit nodes to iterate through
  hitNodes <- xml_children(xml_children(xml_children(xml_children(xmlfile))))
  
  #function to iterate over nodes, check direction, clean and concatenate sequences
  getSequences <- function(x) {
    list <- xml_text(xml_find_all(hitNodes[[x]], xpath = ".//Hsp_hseq"))
    beginning <- as.integer(xml_text(xml_find_all(hitNodes[[x]], xpath = ".//Hsp_hit-from")))
    end <- as.integer(xml_text(xml_find_all(hitNodes[[x]], xpath = ".//Hsp_hit-to")))
    direction <- ifelse(beginning[1]>end[1], "backwards", "normal")
    hit.df <- data.frame(beginning, end, list)
    #ifelse(direction=="normal", hit.df <- arrange(hit.df, beginning, decreasing=FALSE), hit.df <- arrange(hit.df, beginning, decreasing=TRUE))
    if(direction=="normal") {
      hit.df <- arrange(hit.df, beginning, decreasing=FALSE) 
    }  else {
      hit.df <- arrange(hit.df, beginning, decreasing=TRUE)
    } 
    concatenated <- paste(hit.df$list, collapse="")
    gsub("-", "", concatenated)
  }

  
  #############################################
  # workflow
  
  #get all hit ids
  ids <- xml_text(xml_find_all(xmlfile, xpath = "//Hit_def"))
  
  #get all sample names to find out duplicates
  names <- sapply(1:length(ids), splitall) 
  
  #get all sequences
  sequences <- sapply(1:hits, getSequences)
  
  #get all hit lengths; Important: derived from concatenated hits!
  seq_lengths <- nchar(sequences)
  
  #get number of hsps per hit
  hspNum <- sapply(1:hits, getHsps)
  
  #create data frame to search for duplicates etc.
  df <- data.frame(ids, names, seq_lengths, hspNum, sequences)
  #sort data frame by descending sequence length
  df <- df[order(-seq_lengths),]
  
  #find duplicates
  df$duplicate <- duplicated(df$names)
  
  #subset data frame by eliminating duplicates
  dup_rem.df <- df[!(df$duplicate==TRUE),]
  
  
  #############################################
  # create output
  
  #create name for outfile from infile
  sp <- strsplit(url, ".", fixed=TRUE)
  prefix <- sp[[1]][[1]]
  outname <-  paste(prefix, "_contigs.", "fasta", sep="")
  
  #write to fasta
  write.fasta(as.list(dup_rem.df$sequences), dup_rem.df$ids, outname, open="w")  
  
}

sapply(1:length(xml_with_hits), getDatafromXML)
