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
