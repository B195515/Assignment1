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
















