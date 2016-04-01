# fetch-blast-hits
## Description
This R script parses a blast.xml file (or multiple files) and exports the hits into a fasta file. It was created for fetching sequences from an CLC deNovo assembly that match a certain gene locus.

## What it does
The script parses a set of blast-xml-files and extracts all sequences (Hsp_hseq) from each hit and stores them in a fasta file.  
In case there are multiple sequence parts (hsps) per hit, the script checks the orientation and concatenates the hsps together so that they are in the right order. It also deletes all gaps from the blast alignment.
In case there are more that one hit that belong to the same sample, the script checks which of the hits are longer and discardes the shorter ones.

## How to run the script
In order to sort out duplicates, the script needs to recognize which hit belongs to which sample. The sample name must be located after the first "_" in the hit name, like someSortOfID_SAMPLENAME_someOtherInformation.
Here's an example of a potential hit name, with the sample name being "Sebacina-incrustans-SS-16":
C6FNANXX_Sebacina-incrustans-SS-16_AGTGGCAT-TGGCATTC_L002_R1_001_(paired)_trimmed_(paired)_contig_3764

If the hit names are different, the script must be tweaked accordingly in this function:

  #function to get all sample names to find out duplicates
  splitall <- function(x) {
    split.ls <- strsplit(ids[[x]], "_")
    return(split.ls[[1]][[2]])
  }
  
The xml files that need to be parsed must be specified at line 9 of the script. See:
  #get all data
  xml_files <- c("file1.xml", "file2.xml")
