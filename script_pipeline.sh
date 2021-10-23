# set up working space, must initiate code in user's homespace
# all files and subsequent code will be ran while in ~/B195515 directory

cd
echo "Setting up working directory..."
rm -fr B195515
mkdir B195515
cd B195515
mkdir FastaQC
mkdir Refgenome_binary
mkdir Output
mkdir Groups

#!/bin/bash

# Define some variables
AY21=~/B195515/AY21
FASTAQC=~/B195515/FastaQC
FASTQ=~/B195515/AY21/fastq
GENOME=~/B195515/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta
GENOME_ANNOTATED=~/B195515/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed
GENOME_BINARY=~/B195515/Refgenome_binary
OUTPUT=~/B195515/Output
EXPGROUPS=~/B195515/Groups

# Copy original files to B195515 folder
echo
echo "Copying sequence files, reference sequence, and annotated sequence to current folder..."
cp -r /localdisk/data/BPSM/AY21/ .

echo
echo "Performing quality check of sequences..."

# Unzip QC files, output fastqc files to FastaQC
fastqc -q -f fastq -o FastaQC -t 15 $FASTQ/*.fq.gz
gunzip -q $GENOME.gz
unzip -q -o -d $FASTAQC $FASTAQC/\*.zip 

echo
echo "Summarising quality check output for each sequence..."
echo

## Extract header line for SummaryQC
grep -i "Filename\|%GC\|>>" $FASTAQC/100k.C1-1-501_1_fastqc/fastqc_data.txt| grep -v ">>END"| cut -f1| paste -s| sed -E 's/>>//g' > $FASTAQC/SumQC.csv
chmod 755 $FASTAQC/SumQC.csv

## Extract and compile QC summary from each fastaqc file
for i in $FASTAQC/*_fastqc/fastqc_data.txt
do
grep -i "Filename\|%GC\|>>" ${i} | grep -v ">>END" | cut -f2 | paste -s >> $FASTAQC/SumQC.csv
done
cut -f2 $FASTAQC/SumQC.csv > $FASTAQC/SumQC1
cut -f1,3,4,5,6,7,8,9,10,11,12 $FASTAQC/SumQC.csv > $FASTAQC/SumQC2
paste $FASTAQC/SumQC1 $FASTAQC/SumQC2 > $FASTAQC/QCSummary.tsv
cat $FASTAQC/QCSummary.tsv

# Build the binary index files (6 files) all with the specificed prefix
echo
echo "Building binary reference genome sequences..."
bowtie2-build --quiet $GENOME $GENOME_BINARY/refgenome

echo
echo "Successful!"
echo
echo "Now aligning sequences to reference genome with bowtie2..."
echo
echo "This may take a while...so go grab a snack..."
echo

# Align paired-end reads to binary reference genome, then output sorted bam files bypassing sam intermediate, and create matching BAM index file
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bowtie2 -p 15 -x $GENOME_BINARY/refgenome -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz | samtools sort -o ${i}.sorted.bam
done
# Generate index of bam files
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
samtools index -b ${i}.sorted.bam ${i}.sorted.bai
done

echo
echo "Succesful! All bam files are now sorted and indexed." 
echo
echo "Generating bed files and aligning to annotated genome..."

# Converting bam to bed, both are sorted
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bedtools bamtobed -i ${i}.sorted.bam > ${i}.sorted.bed
done
# Align bed files to bed reference, generate sorted count file for each output
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bedtools intersect -a ${i}.sorted.bed -b $GENOME_ANNOTATED -wa -wb | cut -f10 | sort | uniq -c | sed -E 's/^ *//; s/ /\t/' > ${i}.sortedcounts.tsv
done

echo
echo "Now generating average expression levels per gene for each group..."

# Build index file for gene name and description
cut -f4,5 $GENOME_ANNOTATED > genecount.bed
# Join all "sorted counts" file to index file by gene ID
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
join -a 1 -t $'\t' -e '0' -o auto -1 1 -2 2 genecount.bed ${i}.sortedcounts.tsv > output
cat output > genecount.bed
done

# Make a header for counts
echo "GENE NAME" > coltitles
echo "GENE DESCRIPTION" >> coltitles
cat $FASTQ/100k.fqfiles | cut -f1 | grep -v "ID" >> coltitles
# Generate raw counts, compiled
paste -s coltitles > countspergene.tsv
cat genecount.bed >> countspergene.tsv ## move this file

## CLEANUP
mv countspergene.tsv $FASTAQC/QCSummary.tsv $OUTPUT
rm -fr output genecount.bed coltitles $FASTAQC/SumQC1 $FASTAQC/SumQC2 $FASTAQC/SumQC.csv $FASTAQC/SumQC1 $FASTAQC/SumQC2 $FASTAQC/SumQC $FASTAQC/*.zip $FASTAQC/*.zip 

# Sort files by condition
while read ID Sample Replicate Time Treatment End1 End2
do
if [[ $Sample = "Clone1" && $Time = "0" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone1" && $Time = "24" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone1" && $Time = "48" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone1" && $Time = "24" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone1" && $Time = "48" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone2" && $Time = "0" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone2" && $Time = "24" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone2" && $Time = "48" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone2" && $Time = "24" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "Clone2" && $Time = "48" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "WT" && $Time = "0" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "WT" && $Time = "24" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "WT" && $Time = "48" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "WT" && $Time = "24" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
elif [[ $Sample = "WT" && $Time = "48" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
fi
done < $FASTQ/100k.fqfiles

#

echo "All done"

###### THIS IS A TEST AREA #######
