# set up working space, must initiate code in user's homespace
# all files and subsequent code will be ran while in ~/B195515 directory

#-----1. SET UP DIRECTORIES-----#
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

#-----2. SET UP VARIABLES-----#
AY21=~/B195515/AY21
FASTAQC=~/B195515/FastaQC
FASTQ=~/B195515/AY21/fastq
GENOME=~/B195515/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta
GENOME_ANNOTATED=~/B195515/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed
GENOME_BINARY=~/B195515/Refgenome_binary
OUTPUT=~/B195515/Output
EXPGROUPS=~/B195515/Groups

#-----3. COPY ORIGINAL FILES FROM SOURCE TO PWD-----#
echo
echo "Copying sequence files, reference sequence, and annotated sequence to current folder..."
cp -r /localdisk/data/BPSM/AY21/ .

echo
echo "Performing quality check of sequences..."

#-----4. UNZIP REFGENOME, FASTQC; PERFORM FASTQC-----#
fastqc -q -f fastq -o FastaQC -t 15 $FASTQ/*.fq.gz
gunzip -q $GENOME.gz
unzip -q -o -d $FASTAQC $FASTAQC/\*.zip 

echo
echo "Summarising quality check output for each sequence..."
echo

#-----5. SET UP QC-SUMMARY HEADER-----#
## Extract header line for SummaryQC
grep -i "Filename\|%GC\|>>" $FASTAQC/100k.C1-1-501_1_fastqc/fastqc_data.txt| grep -v ">>END"| cut -f1| paste -s| sed -E 's/>>//g' > $FASTAQC/SumQC.csv
chmod 755 $FASTAQC/SumQC.csv

#-----6. COMPILE BRIEF SUMMARY FROM FASTQC OUTPUTS-----#
for i in $FASTAQC/*_fastqc/fastqc_data.txt
do
grep -i "Filename\|%GC\|>>" ${i} | grep -v ">>END" | cut -f2 | paste -s >> $FASTAQC/SumQC.csv
done
cut -f2 $FASTAQC/SumQC.csv > $FASTAQC/SumQC1
cut -f1,3,4,5,6,7,8,9,10,11,12 $FASTAQC/SumQC.csv > $FASTAQC/SumQC2
paste $FASTAQC/SumQC1 $FASTAQC/SumQC2 > $OUTPUT/QCSummary.tsv
cat $OUTPUT/QCSummary.tsv ## move to $OUTPUT

#-----7. BOWTIE2-BUILD: SET UP BINARY INDEX FILES FROM REFGENOME-----#
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

#-----8. BOWTIE2 TO ALIGN PAIRED END READS TO BINARY REFGENOME
#	OUTPUT TO BAM AND CREATE BAM INDEX-----#

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bowtie2 -p 15 -x $GENOME_BINARY/refgenome -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz | samtools sort -o ${i}.sorted.bam
done

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
samtools index -b ${i}.sorted.bam ${i}.sorted.bai
done

echo
echo "Succesful! All bam files are now sorted and indexed." 
echo
echo "Generating bed files and aligning to annotated genome..."

#-----9. BEDTOOLS: BAM -> BED (SORTED) AND ALIGN BED TO ANNOTATED GENOME-----#
#	GENERATE SORTED COUNT FILE FOR EACH ALIGNMENT (2 columns: 1-count, 2-genename)

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bedtools bamtobed -i ${i}.sorted.bam > ${i}.sorted.bed
done

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bedtools intersect -a ${i}.sorted.bed -b $GENOME_ANNOTATED -wa -wb | cut -f10 | sort | uniq -c | sed -E 's/^ *//; s/ /\t/' > ${i}.sortedcounts.tsv
done

echo
echo "Now generating average expression levels per gene for each group..."

#-----10. SORT SAMPLE REPLICATES INTO EXPERIMENTAL GROUP-----#
while read ID Sample Replicate Time Treatment End1 End2
do
if [[ $Sample = "Clone1" && $Time = "0" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "24" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "48" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "24" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "48" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "0" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "24" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "48" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "24" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "48" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "0" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "24" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "48" && $Treatment = "Uninduced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "24" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "48" && $Treatment = "Induced" ]]
then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo $Sample_$Time_$Treatment >> $EXPGROUPS/groupnames
fi
done < $FASTQ/100k.fqfiles
#**code can likely be looped and made flexible for more variables..when enough time..**#


#-----11. COMPILE RAW COUNTS FROM 45 FILES INTO 15 EXPERIMENT GROUPS-----#

ls $EXPGROUPS/*.txt| rev| cut -c 5-| rev| while read groupname
do
cut -f4 $GENOME_ANNOTATED > ${groupname}.groupedcount
done # this sets up count file for each experiment group

for groupname in $(ls $EXPGROUPS/*.txt | rev | cut -c 5-| rev)
do
cat ${groupname}.txt | while read filename
do join -a 1 -t $'\t' -e '0' -o auto -1 1 -2 2 ${groupname}.groupedcount $FASTQ/100k.${filename}.sortedcounts.tsv > output
cat output > ${groupname}.groupedcount
done
done # this reads filenames of each group, then compiles the counts of those filenames into own group (replicates)


#-----12. FIND AVERAGE COUNTS PER GENE PER GROUP, COMPILE GROUPED AVERAGE COUNTS-----#
 
for groupname in $(ls $EXPGROUPS/*.txt | rev | cut -c 5-| rev)
do
awk 'NF > 1 { sum=0; for (i=2; i<=NF; i++) sum+=$i; print $1, (sum/(NF-1)) }' OFS="\t" ${groupname}.groupedcount| sort -k1,1 > ${groupname}.avgcount
done # this finds the average count

cut -f4,5 $GENOME_ANNOTATED > outputavgcount

for groupname in $(ls $EXPGROUPS/*.txt | rev | cut -c 5-| rev)
do
join -a 1 -t $'\t' -1 1 -2 1 outputavgcount ${groupname}.avgcount > output
cat output > outputavgcount
done # this compiles the average count for all experiment groups

echo -e "GENE_NAME\nGENE_DESCRIPTION" > coltitles
cat $EXPGROUPS/groupnames >> coltitles
paste -s coltitles > $OUTPUT/Mean_expression_level.tsv
cat outputavgcount >> $OUTPUT/Mean_expression_level.tsv 










## CLEANUP
rm -fr output coltitles $FASTAQC/SumQC1 $FASTAQC/SumQC2 $FASTAQC/SumQC.csv $FASTAQC/SumQC1 $FASTAQC/SumQC2 $FASTAQC/SumQC $FASTAQC/*.zip $FASTAQC/*.zip

# Sort files by condition
#

echo "All done"

###### THIS IS A TEST AREA #######
