OUTPUT_PREFIX=$1
Threshold=$2

awk -vT=${Threshold} '{if($2==1 && NF==8 && length($3)>=T && length($5)>=T) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_PREFIX.combine.cut.info > $OUTPUT_PREFIX.Combined.linker1.5nd.txt
awk -vT=${Threshold} '{if($2==1 && NF==8 && length($3)>=T && length($5)>=T) print "@"$1"\t"$5"\t""+""\t"$8}' $OUTPUT_PREFIX.combine.cut.info > $OUTPUT_PREFIX.Combined.linker1.3nd.txt
awk -vT=${Threshold} '{if($2==1 && NF==8 && length($3)>=T && length($5)>=T) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_PREFIX.notCombined.1.cut.info  > $OUTPUT_PREFIX.notCombined.linker1.R1.5nd.txt
awk -vT=${Threshold} '{if($2==1 && NF==8 && length($3)>=T && length($5)>=T) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_PREFIX.notCombined.1.cut.info  > $OUTPUT_PREFIX.notCombined.linker1.R1.3nd.txt
awk -vT=${Threshold} '{if($2==1 && NF==8 && length($3)>=T && length($5)>=T) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_PREFIX.notCombined.2.cut.info  > $OUTPUT_PREFIX.notCombined.linker1.R2.5nd.txt
awk -vT=${Threshold} '{if($2==1 && NF==8 && length($3)>=T && length($5)>=T) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_PREFIX.notCombined.2.cut.info  > $OUTPUT_PREFIX.notCombined.linker1.R2.3nd.txt

awk '{if($2==0 && NF==4) print "@"$1"\t"$3"\t""+""\t"$4}' $OUTPUT_PREFIX.notCombined.1.cut.info > $OUTPUT_PREFIX.notCombined.linker0.R1.txt
awk '{if($2==0 && NF==4) print "@"$1"\t"$3"\t""+""\t"$4}' $OUTPUT_PREFIX.notCombined.2.cut.info > $OUTPUT_PREFIX.notCombined.linker0.R2.txt

awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_PREFIX.notCombined.linker0.R2.txt $OUTPUT_PREFIX.notCombined.linker1.R1.5nd.txt > $OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_PREFIX.notCombined.linker0.R2.txt $OUTPUT_PREFIX.notCombined.linker1.R1.3nd.txt > $OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_PREFIX.notCombined.linker1.R2.5nd.txt $OUTPUT_PREFIX.notCombined.linker0.R1.txt > $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_PREFIX.notCombined.linker1.R2.3nd.txt $OUTPUT_PREFIX.notCombined.linker0.R1.txt > $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair.txt

awk '{if(NF==8) print $1"\t"$2"\t"$3"\t"$4}' $OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair.txt > $OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair_R1.txt
awk '{if(NF==8) print $1"\t"$2"\t"$3"\t"$4}' $OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair.txt > $OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair_R1.txt
awk '{if(NF==8) print $5"\t"$6"\t"$7"\t"$8}' $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair.txt > $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair_R2.txt
awk '{if(NF==8) print $5"\t"$6"\t"$7"\t"$8}' $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair.txt > $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair_R2.txt
cat $OUTPUT_PREFIX.Combined.linker1.5nd.txt $OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair_R1.txt $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair_R2.txt > $OUTPUT_PREFIX.total-PET_SE_5nd.txt
cat $OUTPUT_PREFIX.Combined.linker1.3nd.txt $OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair_R1.txt $OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair_R2.txt > $OUTPUT_PREFIX.total-PET_SE_3nd.txt

