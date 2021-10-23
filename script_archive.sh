####### Must run while in the B195515 directory
## To get the header for qc items
grep -i "Filename\|%GC\|>>" ./FastaQC/100k.C1-1-501_1_fastqc/fastqc_data.txt | grep -v ">>END" | cut -f1 | paste -s > ./FastaQC/SumQC.csv
chmod 755 ./FastaQC/SumQC.csv
## To output QC summary from each fastaqc file 
for i in ./FastaQC/*_fastqc/fastqc_data.txt
do
grep -i "Filename\|%GC\|>>" ${i} | grep -v ">>END" | cut -f2 | paste -s >> ./FastaQC/SumQC.csv
done

cut -f2 ./FastaQC/SumQC.csv > ./FastaQC/SumQC1
cut -f1,3,4,5,6,7,8,9,10,11,12 ./FastaQC/SumQC.csv > ./FastaQC/SumQC2
paste ./FastaQC/SumQC1 ./FastaQC/SumQC2 > ./FastaQC/SumQCsorted.csv
cat ./FastaQC/SumQCsorted.csv
# this builds the binary index files (6 files) all with the specelificed prefix
bowtie2-build dir/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta <prefixforoutput> 

# run the alignment for the first file - don't have to gunzip input files
bowtie2 --quiet -x tryp_align1 -1 ./fastq/100k.C1-1-501_1.fq -2 ./fastq/100k.C1-1-501_2.fq \
-S 100k.C1-1-501.sam

##### THIS IS THE FINAL VERSION + LOOPING AND DIRECT OUTPUT TO BAM
for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bowtie2 -p 5 -x trypgenome_ref -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz | samtools view -S -bo ${i}.bam
done


#### making some changes to output sam files into sorted bam

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bowtie2 -p 5 -x trypgenome_ref -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz | samtools sort -S -o ${i}.sorted.bam
done

### editing sam files

bedtools intersect -a 100k.C1-1-501.bed -b fastq/TriTrypDB-46_TcongolenseIL3000_2019.bed -wa -wb > 100k.C1-1-501.wawb.result

### FOR SORTING WAWB FILE
	# get only the column with gene name > sort by gene name > count the number of times gene occurs
bedtools intersect -a 100k.C1-1-501.bed -b fastq/TriTrypDB-46_TcongolenseIL3000_2019.bed -wa -wb | cut -f10 | sort | uniq -c > 100k.C1-1-501.count


## FILE THAT HAS ONLY GENE NAME AND GENE DESCRIPTION: indexgenenamedescription.bed

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bedtools bamtobed -i ${i}.bam > ${i}.bed 
done

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
bedtools intersect -a ${i}.bed -b $AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed -wa -wb | cut -f10 | sort | uniq -c > ${i}.counts.txt
done


## Get index file with Gene description then gene ID 
cut -f4,5 $AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed | awk '{FS="\t"; {print $2,$1}}' > index.bed


######## SHIT WORKS BEAUTIFULLY
cut -f4,5 $AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed > joinindex.bed

for i in $(ls $FASTQ/*.fq.gz | rev | cut -c 9- | rev | uniq)
do
join -a 1 -e -o -1 1 -2 2 joinindex.bed ${i}.sortedcounts.tsv > output
cat output > joinindex.bed
done

###### JOIN AT FIELD1 FILE1 and FIELD2 FILE2, show unpairabble lines from file1, empty, obey format
join -a 1 -e -o -1 1 -2 2 joinindex2.tsv 100k.WT-6-555_1.sortedcounts.txt > joined.txt

join -a 1 -e -o -1 1 -2 2 joined.txt 

### sample code
while read line; do
    join -a 1 file1 "$line" > output
    cat output > file1
    done < list

#### get the files in groups

while read ID Sample Replicate Time Treatment End1 End2
do
elif [[ $Sample == "Clone1" ]]  &&  [[ $Time == "0" ]] && [[ $Treatment == "Uninduced" ]]
then echo -e "$filename" >> Clone1_t0_uninduced.txt
elif [[ $Sample == "Clone1"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> Clone1_t24_uninduced.txt
elif [[ $Sample == "Clone1"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> Clone1_t48_uninduced.txt
elif [[ $Sample == "Clone1"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then echo -e "$filename" >> Clone1_t24_induced.txt
elif [[ $Sample == "Clone1"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then echo -e "$filename" >> Clone1_t48_induced.txt
elif [[ $Sample == "Clone2"]] && [[$Time == "0"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> Clone2_t0_uninduced.txt
elif [[ $Sample == "Clone2"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> Clone2_t24_uninduced.txt
elif [[ $Sample == "Clone2"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> Clone2_t48_uninduced.txt
elif [[ $Sample == "Clone2"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then echo -e "$filename" >> Clone2_t24_induced.txt
elif [[ $Sample == "Clone2"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then echo -e "$filename" >> Clone2_t48_induced.txt
elif [[ $Sample == "WT"]] && [[$Time == "0"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> WT_t0_uninduced.txt
elif [[ $Sample == "WT"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> WT_t24_uninduced.txt
elif [[ $Sample == "WT"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then echo -e "$filename" >> WT_t48_uninduced.txt
elif [[ $Sample == "WT"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then echo -e "$filename" >> WT_t24_induced.txt
elif [[ $Sample == "WT"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then echo -e "$filename" >> WT_t48_induced.txt
fi
done < $FASTQ/100k.fqfiles

##########################
while read ID Sample Replicate Time Treatment End1 End2
do
scorebin=1
if [[ $Sample == "Clone1" ]] && [[ $Time == "0" ]] && [[ $Treatment == "Uninduced" ]]
then scorebin=1
elif [[ $Sample == "Clone1"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then scorebin=2
elif [[ $Sample == "Clone1"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then scorebin=3
elif [[ $Sample == "Clone1"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then scorebin=4
elif [[ $Sample == "Clone1"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then scorebin=5
elif [[ $Sample == "Clone2"]] && [[$Time == "0"]] && [[$Treatment == "Uninduced" ]]
then scorebin=6
elif [[ $Sample == "Clone2"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then scorebin=7
elif [[ $Sample == "Clone2"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then scorebin=8
elif [[ $Sample == "Clone2"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then scorebin=9
elif [[ $Sample == "Clone2"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then scorebin=10
elif [[ $Sample == "WT"]] && [[$Time == "0"]] && [[$Treatment == "Uninduced" ]]
then scorebin=11
elif [[ $Sample == "WT"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then scorebin=12
elif [[ $Sample == "WT"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then scorebin=13
elif [[ $Sample == "WT"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then scorebin=14
elif [[ $Sample == "WT"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then scorebin=15
fi
filename=$(echo -e "$ID")
case $scorebin in
1)
echo -e "$filename" >> Clone1_t0_uninduced.txt ;;
2)
echo -e "$filename" >> Clone1_t24_uninduced.txt ;;
3)
echo -e "$filename" >> Clone1_t48_uninduced.txt ;;
4)
echo -e "$filename" >> Clone1_t24_induced.txt ;;
5)
echo -e "$filename" >> Clone1_t48_induced.txt ;;
6)
echo -e "$filename" >> Clone2_t0_uninduced.txt ;;
7)
echo -e "$filename" >> Clone2_t24_uninduced.txt ;;
8)
echo -e "$filename" >> Clone2_t48_uninduced.txt ;;
9)
echo -e "$filename" >> Clone2_t24_induced.txt ;;
10)
echo -e "$filename" >> Clone2_t48_induced.txt ;;
11)
echo -e "$filename" >> WT_t0_uninduced.txt ;;
12)
echo -e "$filename" >> WT_t24_uninduced.txt ;;
13)
echo -e "$filename" >> WT_t48_uninduced.txt ;;
14)
echo -e "$filename" >> WT_t24_induced.txt ;;
15)
echo -e "$filename" >> WT_t48_induced.txt ;;
esac
done < $FASTQ/100k.fqfiles
##################################

##################################
while read ID Sample Replicate Time Treatment End1 End2
do
if [[ $Sample == "Clone1" ]] && [[ $Time == "0" ]] && [[ $Treatment == "Uninduced" ]]
then scorebin=1
elif [[ $Sample == "Clone1"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then scorebin=2
elif [[ $Sample == "Clone1"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then scorebin=3
elif [[ $Sample == "Clone1"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then scorebin=4
elif [[ $Sample == "Clone1"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then scorebin=5
elif [[ $Sample == "Clone2"]] && [[$Time == "0"]] && [[$Treatment == "Uninduced" ]]
then scorebin=6
elif [[ $Sample == "Clone2"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then scorebin=7
elif [[ $Sample == "Clone2"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then scorebin=8
elif [[ $Sample == "Clone2"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then scorebin=9
elif [[ $Sample == "Clone2"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then scorebin=10
elif [[ $Sample == "WT"]] && [[$Time == "0"]] && [[$Treatment == "Uninduced" ]]
then scorebin=11
elif [[ $Sample == "WT"]] && [[$Time == "24"]] && [[$Treatment == "Uninduced" ]]
then scorebin=12
elif [[ $Sample == "WT"]] && [[$Time == "48"]] && [[$Treatment == "Uninduced" ]]
then scorebin=13
elif [[ $Sample == "WT"]] && [[$Time == "24"]] && [[$Treatment == "Induced" ]]
then scorebin=14
elif [[ $Sample == "WT"]] && [[$Time == "48"]] && [[$Treatment == "Induced" ]]
then scorebin=15
fi

