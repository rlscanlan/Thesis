#https://cran.r-project.org/web/packages/fastqcr/fastqcr.pdf

#install.packages("devtools")
#install.packages("rvest")
library(devtools)
library(fastqcr)
library(dplyr)

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/fastqcr")

fastqc_install(dest.dir = "~/bin")
#wd for mac
wd = setwd("/Volumes/Bex5TB/Sen23")
#wd for linux
wd = setwd("/media/c0068011/Bex5TB/Sen23")

wd = setwd("/Volumes/Bex5TB/Sen23/GSE179465/")
directories = list.dirs(full.names = FALSE, recursive = FALSE)
print(directories)


for(i in directories){
  print(paste0(i, "fastq"))
  fastqc(
    fq.dir = paste0(i, "/"),
    qc.dir = paste0(i, "/fastqc"),
    threads =4,
    fastqc.path = ("~/bin/FastQC/fastqc")
  )
}
#qc creates aggregate of ALL samples within directory so all GSEs
qc = qc_aggregate(qc.dir = ".", progressbar = TRUE)
summary(qc)



#to get .html multiQC report, move over to terminal and make sure in python
#run: pip install multiqc
#setwd to location. On mac: cd /Volumes/Bex5TB/Sen23
#run: multiqc .


####cutadapt####
#https://cutadapt.readthedocs.io/en/v4.4/guide.html
#cutadapt in terminal:
# conda activate cutadaptenv
# to check it worked run: cutadapt --version
#trim based on quality. set quality at 25 as anything above a phred score of 20 is often accepted as acceptable
#as this means there is only a 1% chance of error, therefore this is reduced at a phred of 25
#Example run which takes the file (second named file) and creates the output file (first named file) with a q of 25:
#cutadapt -q 25 -o SRR16848303_GSM5678364_IMR90_ER_RAS_d5_rep3_Homo_sapiens_RNA-Seq_q25.fastq.gz SRR16848303_GSM5678364_IMR90_ER_RAS_d5_rep3_Homo_sapiens_RNA-Seq.fastq.gz
#can run more than one line at once
