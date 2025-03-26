#in terminal: conda activate salmon

#create index. if salmon is updated, a new index will need creating in new salmon version
#first get the fasta files from following link https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/

#use wget to get if know link. On October 4th 2023 the links were as below
#wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
#wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

#zcat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz >Homo_sapiens_h38_110_cdna_ncrna.fa

#in terminal run the following line to create index
#salmon index -t Homo_sapiens_h38_110_cdna_ncrna.fa -i Homo_sapiens_h38_110_cdna_ncrna_quasi_index_Salmon1.10.1


#now move onto running samples through salmon
#set wd. in my case here I need: cd /media/c0068011/Bex5TB/Sen23
#in terminal: bash quant_salmon_Sen23samples.sh

#set directory to Sen23
index="/media/c0068011/Bex5TB/Salmon/Homo_sapiens_h38_110_cdna_ncrna_quasi_index_Salmon1.10.1"

####for paired library layout####
for fn in ./Downloaded_on_Oct23/GSE175533/SRR14646{293..322};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -1 ${fn}/${samp}*_1.fastq.gz \
         -2 ${fn}/${samp}*_2.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE175533/quants/${samp}_quant
done 


for fn in ./Downloaded_on_Oct23/GSE175533/SRR14646{353..370};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -1 ${fn}/${samp}*_1.fastq.gz \
         -2 ${fn}/${samp}*_2.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE175533/quants/${samp}_quant
done 


for fn in ./Downloaded_on_Oct23/GSE196610/SRR179887{12..35};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -1 ${fn}/${samp}*_1.fastq.gz \
         -2 ${fn}/${samp}*_2.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE196610/quants/${samp}_quant
done 


for fn in ./Downloaded_on_Oct23/GSE221104/SRR22756{192..200};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -1 ${fn}/${samp}*_1.fastq.gz \
         -2 ${fn}/${samp}*_2.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE221104/quants/${samp}_quant
done 


for fn in ./Downloaded_on_Oct23/GSE225095/SRR2341816{0..5};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -1 ${fn}/${samp}*_1.fastq.gz \
         -2 ${fn}/${samp}*_2.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE225095/quants/${samp}_quant
done


for fn in ./Downloaded_on_Oct23/GSE240226/SRR255574{88..99};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -1 ${fn}/${samp}*_1.fastq.gz \
         -2 ${fn}/${samp}*_2.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE240226/quants/${samp}_quant
done


####single-end library layout####
for fn in ./Downloaded_on_Oct23/GSE224070/SRR23272{462..503};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -r ${fn}/${samp}*.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE224070/quants/${samp}_quant
done

for fn in ./Downloaded_on_Oct23/GSE224071/SRR232724{44..61};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -r ${fn}/${samp}*.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./SingDownloaded_on_Oct23le/GSE224071/quants/${samp}_quant
done

for fn in ./Downloaded_on_Oct23/GSE234417/SRR248720{01..16};
do
samp=`basename ${fn}*`
echo "Processing sample ${samp}"
salmon quant -i ${index} -l A \
         -r ${fn}/${samp}*.fastq.gz \
         -p 8 --gcBias --seqBias --validateMappings -o ./Downloaded_on_Oct23/GSE234417/quants/${samp}_quant
done


#for fn in ./Single/GSE179465/Trimmed/SRR150437{68..73};
#do
#samp=`basename ${fn}*`
#echo "Processing sample ${samp}"
#salmon quant -i ${index} -l A \
#         -r ${fn}/${samp}*.fastq.gz \
#         -p 8 --gcBias --seqBias --validateMappings -o ./Single/GSE179465/Trimmed/quants/${samp}_quant
#done


####technical replicates in single-end study####
#gse155903
#first need to combine the fastq files from the technical repeats
#combined first because of this post https://support.bioconductor.org/p/9153803/#9153994
#I've done this using the following file:
#/media/c0068011/Bex5TB/Sen23/Single/GSE155903/combine_tech_repeats_GSE155903.sh (run in terminal using bash then file name)
#then run salmon like normal for the library layout design


#for fn in ./Single/GSE155903/Combined/GSM47154{45..56};
#do
#samp=`basename ${fn}*`
#echo "Processing sample ${samp}"
#salmon quant -i ${index} -l A \
#         -r ${fn}/${samp}*.fastq.gz \
#         -p 8 --gcBias --seqBias --validateMappings -o ./Single/GSE155903/quants/${samp}_quant
#done