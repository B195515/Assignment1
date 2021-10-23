## COUNT MEAN PER ROW FOR COLUMN 3 ONWARDS; PRINT OUT COL1 AND COL2 AND AVERAGE


### Script 001
awk 'NR == 1 {print $1, $2, "Average"; next } ## f1 is gene name, f2 is gene description
NF > 2 { sum=0; for (i=3; i<=NF; i++) sum+=$i; print $1, $2, (sum/(NF-2)) }' OFS="\t" countmean.xls > outputfile
	### Script 001 output ## does not give the correct mean!
GENE    NAME    Average
TcIL3000.A.H_000005000.1        Retrotransposon 1.42857
TcIL3000.A.H_000005100.1        C-5     3.83333
TcIL3000.A.H_000005200.1        hypothetical    0


### Script 002
awk 'NR == 1 {FS="\t"; print $1, "Average"; next }
NF > 2 { sum=0; for (i=3; i<=NF; i++) sum+=$i; print $1, (sum/(NF-2)) }' OFS="\t" countmean.xls > outputfile
	### Script 002 output - gives correct mean, but column header looks weird..	
GENE    Average
TcIL3000.A.H_000005000.1        3.33333
TcIL3000.A.H_000005100.1        7.66667
TcIL3000.A.H_000005300.1        9.66667

### Script 003
awk 'NR == 1 {FS="\t"; print "GENE_NAME", "Average"; next }
NF > 2 { sum=0; for (i=3; i<=NF; i++) sum+=$i; print $1, (sum/(NF-2)) }' OFS="\t" test_countmean > test_output
	### Script 003 output - correct mean, forced column header
	GENE_NAME       Average
	TcIL3000.A.H_000005000.1        3.33333
	TcIL3000.A.H_000005100.1        7.66667
	TcIL3000.A.H_000005200.1        0
	TcIL3000.A.H_000005300.1        9.66667

################
while read filename
do
if $filename 
awk 'NR == 1 {print $1, $2, "Average"; next }
NF > 2 { sum=0; for (i=3; i<=NF; i++) sum+=$i; print $1, $2, (sum/(NF-2)) }' countmean.xls
	#### need to put some code here
done < $EXPGROUPS/*.txt
