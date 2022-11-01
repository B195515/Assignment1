#-----set up working space, must initiate script in user's homespace-----#
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

#-----2. SET UP VARIABLES PATHS-----#
AY21=~/B195515/AY21
FASTAQC=~/B195515/FastaQC
FASTQ=~/B195515/AY21/fastq
GENOME=~/B195515/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta
GENOME_ANNOTATED=~/B195515/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed
GENOME_BINARY=~/B195515/Refgenome_binary
OUTPUT=~/B195515/Output
EXPGROUPS=~/B195515/Groups

#-----3. COPY ORIGINAL FILES FROM SOURCE TO PWD-----#
echo "Copying sequence files, reference sequence, and annotated sequence to current folder..."
cp -r /localdisk/data/BPSM/AY21/ .

#-----4. QUALITY CHECKING-----#
echo "Performing quality check of sequences..."

fastqc -q -f fastq -o FastaQC -t 15 $FASTQ/*.fq.gz
gunzip -q $GENOME.gz
unzip -q -o -d $FASTAQC $FASTAQC/\*.zip 
echo "Summarising quality check output for each sequence..."

	# this generates QC summary file header line
grep -i "Filename\|%GC\|>>" $FASTAQC/100k.C1-1-501_1_fastqc/fastqc_data.txt| grep -v ">>END"| cut -f1| paste -s| sed -E 's/>>//g' > $FASTAQC/SumQC.csv
chmod 755 $FASTAQC/SumQC.csv 

	# this generates a summary of QC from fastqc program
for i in $FASTAQC/*_fastqc/fastqc_data.txt; do
grep -i "Filename\|%GC\|>>" ${i} | grep -v ">>END" | cut -f2 | paste -s >> $FASTAQC/SumQC.csv; done
cut -f2 $FASTAQC/SumQC.csv > $FASTAQC/SumQC1
cut -f1,3,4,5,6,7,8,9,10,11,12 $FASTAQC/SumQC.csv > $FASTAQC/SumQC2
paste $FASTAQC/SumQC1 $FASTAQC/SumQC2 > $OUTPUT/QCSummary.tsv

#-----5. ALIGNING SEQUENCES TO REFERENCE GENOME-----#

	# this generates binary files from ref genome fasta file, outputs files with the same suffix specified
echo "Building binary reference genome sequences..."
bowtie2-build --quiet $GENOME $GENOME_BINARY/refgenome

echo -e "Successful!\nNow aligning sequences to reference genome with bowtie2..."

	# this outputs alignment from each sequence in the sample that aligned to reference genome, producing sorted bam files
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq); do
bowtie2 -p 15 -x $GENOME_BINARY/refgenome -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz | samtools sort -o ${i}.sorted.bam; done

	# this outputs an index for the bam
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq); do
samtools index -b ${i}.sorted.bam ${i}.sorted.bai; done 

echo -e "Succesful! All bam files are now sorted and indexed.\nGenerating bed files and aligning to annotated genome..."

#-----6. INDEX ALIGNED SEQUENCES TO ANNOTATED REFERENCE GENOME BY BEDTOOLS-----#
	# this converts bam to bed (presorted)
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq); do
bedtools bamtobed -i ${i}.sorted.bam > ${i}.sorted.bed; done
	
	# this outputs regions in B that aligned with A, for each query in A, then produces a sorted count file for each sample 
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq); do
bedtools intersect -a ${i}.sorted.bed -b $GENOME_ANNOTATED -wa -wb | cut -f10 | sort | uniq -c | sed -E 's/^ *//; s/ /\t/' > ${i}.sortedcounts.tsv; done
	#-----MODIFY/COMBINE BEDTOOLS PARAMETERS FOR OTHER TYPES OF ANALYSES-----#

	# this generates compilation of raw expression level for each sample (not grouped)
	# empty fields are output as '0' rather than 'NA' to simplify calculations with awk later...
cut -f4,5 $GENOME_ANNOTATED > outputrawexpression
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq); do
join -a 1 -t $'\t' -e '0' -o auto -1 1 -2 2 outputrawexpression ${i}.sortedcounts.tsv > output
cat output > outputrawexpression; done 

	# this compiles the raw expression final output file
echo -e "GENE_NAME\nGENE_DESCRIPTION" > coltitles
cat $FASTQ/100k.fqfiles | cut -f1 | grep -v "ID" >> coltitles
paste -s coltitles > $OUTPUT/Raw_expression_level.tsv
cat outputrawexpression >> $OUTPUT/Raw_expression_level.tsv 

echo -e "Successful!\nNow generating average expression levels per gene for each group..."

#-----7. SORT SAMPLE REPLICATES INTO EXPERIMENTAL GROUP-----#

	# this sorts sample replicates into experiment groups that can be referred to later
	# NOTE: this section can likely be looped and made flexible for more variables in future.. #
while read ID Sample Replicate Time Treatment End1 End2; do
if [[ $Sample = "Clone1" && $Time = "0" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "24" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "48" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "24" && $Treatment = "Induced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone1" && $Time = "48" && $Treatment = "Induced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "0" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "24" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "48" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "24" && $Treatment = "Induced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "Clone2" && $Time = "48" && $Treatment = "Induced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "0" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "24" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "48" && $Treatment = "Uninduced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "24" && $Treatment = "Induced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames
elif [[ $Sample = "WT" && $Time = "48" && $Treatment = "Induced" ]]; then
echo $ID >> $EXPGROUPS/${Sample}_${Time}_${Treatment}.txt
echo ${Sample}_${Time}_${Treatment} >> $EXPGROUPS/groupnames; fi; done < $FASTQ/100k.fqfiles 
echo "Sorting group replicates and generating average expression levels..."

	# this sets up count file for each experiment group
ls $EXPGROUPS/*.txt| rev| cut -c 5-| rev| while read groupname; do
cut -f4 $GENOME_ANNOTATED > ${groupname}.groupedcount; done 

	# empty fields are output as '0' instead of 'NA' to simplify downstream awk script
	# this reads filenames of each group, then compiles the counts of those filenames into own group (replicates)
for groupname in $(ls $EXPGROUPS/*.txt | rev | cut -c 5-| rev); do
cat ${groupname}.txt | while read filename; do
join -a 1 -t $'\t' -e '0' -o auto -1 1 -2 2 ${groupname}.groupedcount $FASTQ/100k.${filename}.sortedcounts.tsv > output
cat output > ${groupname}.groupedcount; done; done 

#-----8. GENERATE AVERAGE COUNTS PER GENE-----#

	# this finds the average count
for groupname in $(ls $EXPGROUPS/*.txt | rev | cut -c 5-| rev); do
awk 'NF > 1 { sum=0; for (i=2; i<=NF; i++) sum+=$i; print $1, (sum/(NF-1)) }' OFS="\t" ${groupname}.groupedcount| sort -k1,1 > ${groupname}.avgcount; done 

	# this compiles the average count for all experiment groups
cut -f4,5 $GENOME_ANNOTATED > outputavgcount
for groupname in $(ls $EXPGROUPS/*.txt | rev | cut -c 5-| rev); do
join -a 1 -t $'\t' -1 1 -2 1 outputavgcount ${groupname}.avgcount > output
cat output > outputavgcount; done 

	# this sets up table to view mean expression level per group (n=3/group) 
echo -e "GENE_NAME\nGENE_DESCRIPTION" > coltitles
cat $EXPGROUPS/groupnames| sort| uniq >> coltitles
paste -s coltitles > $OUTPUT/Mean_expression_level.tsv
cat outputavgcount >> $OUTPUT/Mean_expression_level.tsv 

#-----9. GENE EXPRESSION FOLD CHANGE OVER CONTROL - VARIOUS GROUP COMPARISONS-----#
echo "Generating fold change expression for various group comparison..."

	# this groups count data by time=48hours
awk 'NR==1 {
for (i=1;i<=NF;++i) if ($i=="WT_48_Uninduced") { n=i }
} {for (i=1;i<=NF;++i) if ($i=="WT_48_Induced") { p=i }
} {for (i=1;i<=NF;++i) if ($i=="Clone1_48_Uninduced") { q=i } 
} {for (i=1;i<=NF;++i) if ($i=="Clone1_48_Induced") { r=i } 
} {for (i=1;i<=NF;++i) if ($i=="Clone2_48_Uninduced") { s=i } 
} {for (i=1;i<=NF;++i) if ($i=="Clone2_48_Induced") { v=i; break }} {FS=IFS="\t"; print $1, $2, $n, $p, $q, $r, $s, $v }' OFS="\t" $OUTPUT/Mean_expression_level.tsv > output

	# this calculates fold change using WT_48_Uninduced as baseline of "1"
awk 'NR==1 {
FS="\t"; print "GENE_NAME", "GENE_DESCRIPTION", "FOLDCHANGE_WT_48Ind_vs_WT_48Unind", "FOLDCHANGE_C1_48Unind_vs_WT_48Unind", "FOLDCHANGE_C1_48Ind_vs_WT_48Unind", "FOLDCHANGE_C2_48Unind_vs_WT_48Unind", "FOLDCHANGE_C2_48Ind_vs_WT_48Unind"; next} NR>1 {
if ($3>0) {
FS="\t"; $8=($4/$3); $9=($5/$3); $10=($6/$3); $11=($7/$3); $12=($8/$3); print $1, $2, $8, $9, $10, $11, $12} else {
print $1, $2, "NA", "NA", "NA", "NA", "NA"}
}' OFS="\t" output > $OUTPUT/Foldchange_48h.tsv

#----- ADD OTHER GROUP-WISE COMPARISONS HERE-----#

#-----CLEANUP-----#
echo "Performing cleanup of temp files..."
rm -fr output outputavgcount outputrawexpression coltitles $FASTAQC/SumQC1 $FASTAQC/SumQC2 $FASTAQC/SumQC.csv $FASTAQC/*.zip
echo "Successful!"

#-----END NOTES-----#
SEQTECHNIQUE=RNAseq
EXPORGANISM="Trypanosoma congolense"
EXPTECHNIQUE="synthetic RNAi construct"
OUTPUTDIRNAME="B195515/Output"
echo -e "### END MESSAGE ###\nBasic analysis for paired-end sequence reads generated from $SEQTECHNIQUE on $EXPORGANISM manipulation by $EXPTECHNIQUE is complete.\nFinal result files are located in $OUTPUTDIRNAME and contains tab-delimited files summarising quality of sequence output, raw data counts, averaged gene expression level, and group-wise comparisons.\nFollowing are the tools used for this analysis:"
echo
gunzip --version
echo
fastqc --version
echo
bowtie2 --version
echo
bedtools --version
echo
samtools --version
echo
echo "End"
