
## python3 environment needed !!
## Usage: Usage: sh Run.sh ${prefix} linker_file Int[length(linker)*0.9]
################################################################################################################################################################################

OUTPUT_PREFIX=$1
Linker_file=$2

OUTPUT_DIRECTORY1='Step1_extract_linker_record'
OUTPUT_DIRECTORY2='Step2_pre_process'
OUTPUT_DIRECTORY3='Step3_Align'
OUTPUT_DIRECTORY4='Step4_recognize_real_chimeric'
OUTPUT_DIRECTORY5='Step5_Intra-Molecular_Cluster'
OUTPUT_DIRECTORY6='Step6_Inter-Molecular_Cluster'

mkdir $OUTPUT_DIRECTORY1
mkdir $OUTPUT_DIRECTORY2
mkdir $OUTPUT_DIRECTORY3
mkdir $OUTPUT_DIRECTORY4
mkdir $OUTPUT_DIRECTORY5
mkdir -p $OUTPUT_DIRECTORY6/sim

################################################################################################################################################################################

IN1=${OUTPUT_PREFIX}_R1.fq
IN2=${OUTPUT_PREFIX}_R2.fq

Nsim=100
THREADS='8'
HISAT2_GENOME_INDEX=""
genome_info="genome.element.bed"

gunzip ./*gz
module load BEDTools/2.27

################################################################################################################################################################################



################################################################################################
################################## 1. linker identification ####################################
################################################################################################

####################################### 1.0 trim ###############################################
# java -jar \
# ~/Trimmomatic-0.39/trimmomatic-0.39.jar \
# PE \
# -phred33 \
# -threads ${THREADS} \
# ${IN1} ${IN2} \
# ${OUT}_pair_R1.fq \
# ${OUT}_unpair_R1.fq \
# ${OUT}_pair_R2.fq \
# ${OUT}_unpair_R2.fq \
# ILLUMINACLIP:${Trimmomatichome}/adapters/TruSeq3-PE.fa:2:30:7:8:true \
# LEADING:10 \
# TRAILING:10 \
# SLIDINGWINDOW:4:15 \
# MINLEN:30

##################################### 1.1 flash Combine #############################################
flash --threads $THREADS -M 145 -O --output-prefix $OUTPUT_PREFIX.flash --output-directory $OUTPUT_DIRECTORY1 ${IN1} ${IN2}


##################################### 1.2 linker recognization ######################################
cutadapt -b file:$Linker_file -n 8 -e 0.2 --no-indels -o $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.noLinker       --info-file $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.Linker_info       --discard -O 34 $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.flash.extendedFrags.fastq > $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.stat
cutadapt -b file:$Linker_file -n 8 -e 0.2 --no-indels -o $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.noLinker --info-file $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.Linker_info --discard -O 34 $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.flash.notCombined_1.fastq > $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.stat
cutadapt -b file:$Linker_file -n 8 -e 0.2 --no-indels -o $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.noLinker --info-file $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.Linker_info --discard -O 34 $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.flash.notCombined_2.fastq > $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.stat


##################################### 1.3 linker info sta ###########################################
perl ${scriptdir}/linker_info_sta.pl $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.Linker_info       $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.cut.info       > $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.cut.info.log
perl ${scriptdir}/linker_info_sta.pl $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.Linker_info $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.cut.info > $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.cut.info.log
perl ${scriptdir}/linker_info_sta.pl $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.Linker_info $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.cut.info > $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.cut.info.log


##################################### 1.4 extract valid data ########################################
awk '{if($2==1 && NF==8 && length($3)>=13 && length($5)>=13) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.cut.info > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.Combined.linker1.5nd.txt
awk '{if($2==1 && NF==8 && length($3)>=13 && length($5)>=13) print "@"$1"\t"$5"\t""+""\t"$8}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.combine.cut.info > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.Combined.linker1.3nd.txt
awk '{if($2==1 && NF==8 && length($3)>=13 && length($5)>=13) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.cut.info  > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R1.5nd.txt
awk '{if($2==1 && NF==8 && length($3)>=13 && length($5)>=13) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.cut.info  > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R1.3nd.txt
awk '{if($2==1 && NF==8 && length($3)>=13 && length($5)>=13) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.cut.info  > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R2.5nd.txt
awk '{if($2==1 && NF==8 && length($3)>=13 && length($5)>=13) print "@"$1"\t"$3"\t""+""\t"$6}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.cut.info  > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R2.3nd.txt

awk '{if($2==0 && NF==4) print "@"$1"\t"$3"\t""+""\t"$4}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.1.cut.info > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0.R1.txt
awk '{if($2==0 && NF==4) print "@"$1"\t"$3"\t""+""\t"$4}' $OUTPUT_DIRECTORY1/$OUTPUT_PREFIX.notCombined.2.cut.info > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0.R2.txt

awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0.R2.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R1.5nd.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0.R2.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R1.3nd.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R2.5nd.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0.R1.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair.txt
awk 'NR==FNR{a[$1]=$0}NR>FNR{print $0"\t"a[$1]}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1.R2.3nd.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0.R1.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair.txt

awk '{if(NF==8) print $1"\t"$2"\t"$3"\t"$4}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair_R1.txt
awk '{if(NF==8) print $1"\t"$2"\t"$3"\t"$4}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair_R1.txt
awk '{if(NF==8) print $5"\t"$6"\t"$7"\t"$8}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair_R2.txt
awk '{if(NF==8) print $5"\t"$6"\t"$7"\t"$8}' $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair_R2.txt
cat $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.Combined.linker1.5nd.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-5nd.linker0-R2.pair_R1.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-5nd.pair_R2.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_5nd.txt
cat $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.Combined.linker1.3nd.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker1-R1-3nd.linker0-R2.pair_R1.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.notCombined.linker0-R1.linker1-R2-3nd.pair_R2.txt > $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_3nd.txt


##################################### 1.5  Reverse complementation of the 3' valid sequence #########
sed -i "s/\t/\n/g" $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_3nd.txt
python ${scriptdir}/sequence_RC.py $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_3nd.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_3nd_FXHBcorrected.fq

cp $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_5nd.txt $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_5nd.fq
sed -i "s/\t/\n/g" $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_5nd.fq

gzip $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_5nd.fq
gzip $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_3nd_FXHBcorrected.fq


##################################### 1.6  Duplication remove #######################################
python ${scriptdir}/remove_duplication.py \
-r1 $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_5nd.fq.gz \
-r2 $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.total-PET_SE_3nd_FXHBcorrected.fq.gz \
-o1 $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.part1.fq \
-o2 $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.part2.fq >> $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.rmDup.log





################################################################################################
################################### 2. sequence alignment ######################################
################################################################################################

##################################### 2.1 Pair-end RNA sequences alignment ##########################
hisat2 -p ${THREADS} -x ${HISAT_GENOME_INDEX} -U $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.part1.fq -S $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1.sam --summary-file $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_hisat2_aln_out.txt
hisat2 -p ${THREADS} -x ${HISAT_GENOME_INDEX} -U $OUTPUT_DIRECTORY2/$OUTPUT_PREFIX.part2.fq -S $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2.sam --summary-file $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_hisat2_aln_out.txt


##################################### 2.2 extract Pair-end valid alignment record ###################
cat $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1.sam | grep 'NH:i:1$' > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_unique.sam
cat $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2.sam | grep 'NH:i:1$' > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_unique.sam
awk '{if($5>=20) print $0}' $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_unique.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_unique.q20.sam
awk '{if($5>=20) print $0}' $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_unique.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_unique.q20.sam

awk '{print $1}' $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_unique.q20.sam >> $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.aln_list.txt
awk '{print $1}' $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_unique.q20.sam >> $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.aln_list.txt
sort $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.aln_list.txt |uniq -c | sort -nrk 1 > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.aln_list_sta.txt                 # 5,998,627
awk '{if($1==2) print $2}' $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.aln_list_sta.txt > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.pair_aln_uniqueID_list.txt    # 1,324,583

grep -Fwf $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.pair_aln_uniqueID_list.txt $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_unique.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_PET_unique_aln.sam
grep -Fwf $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.pair_aln_uniqueID_list.txt $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_unique.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_PET_unique_aln.sam


##################################### 2.3 Valid alignment data combination ##########################
samtools view -H $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.header.sam

cat $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part1_PET_unique_aln.sam $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part2_PET_unique_aln.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part_PET_unique_aln.sam
cat $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.header.sam $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.part_PET_unique_aln.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.PET_unique_aln.sam
samtools view -bS $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.PET_unique_aln.sam > $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.PET_unique_aln.bam
samtools sort -n -o $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.PET_unique_aln_sort.bam $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.PET_unique_aln.bam





#####################################################################################################
################################### 3. Cluster Identification #######################################
#####################################################################################################

Mscriptdir=/public/home/spluan/other/Scripts/02_2023_DIY/RICseq/RICATK_lab/modified_4_HiR
#################################### 3.1 Inter- and Intra-molecular interacion classification #######
python ${scriptdir}/class_bam_PairTag.py \
-g ${genome_info} \
-i $OUTPUT_DIRECTORY3/$OUTPUT_PREFIX.PET_unique_aln_sort.bam \
-o $OUTPUT_DIRECTORY4/$OUTPUT_PREFIX.genome > $OUTPUT_DIRECTORY4/$OUTPUT_PREFIX.genome.class.log


######################################### 3.2.1 convert sam/bam into bedpe file #####################
python ${scriptdir}/IntraSam2bedpe.py -i $OUTPUT_DIRECTORY4/$OUTPUT_PREFIX.genome.interaction.intraMolecular.bam -o $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.bedpe
LANG=C sort --parallel=${THREADS} -k1,1 -k8,8 -k2,2n -k3,3n -k5,5n -k6,6 $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.bedpe > $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.sorted.bedpe


######################################### 3.2.2 cluster calling #####################################
java -jar ${scriptdir}/ClusterCalling.jar call $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.sorted.bedpe $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.cluster


######################################### 3.2.3 call high score interacion ##########################
java -jar ${scriptdir}/ClusterCalling.jar count $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.sorted.bedpe $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.cluster $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.cluster.count
Rscript ${scriptdir}/Hypergeometric.r $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.cluster.count $OUTPUT_DIRECTORY5/$OUTPUT_PREFIX.genome.interaction.intraMolecular.cluster.count.pvalue



#################################### 3.3.1 convert sam/bam into bedpe file ##########################
python ${Mscriptdir}/Interbam2bedpe.py -i $OUTPUT_DIRECTORY4/$OUTPUT_PREFIX.genome.interaction.interMolecular.bam -o $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.genome.interaction.interMolecular.bedpe
pigz -p 8 $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.genome.interaction.interMolecular.bedpe


#################################### 3.3.2 get inter-interaction ####################################
python ${scriptdir}/getInteraction.py \
-i $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.genome.interaction.interMolecular.bedpe.gz \
-o $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.total.bedpe.interaction.txt \
-l ${genome_info}

pigz -p ${THREADS} $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.total.bedpe.interaction.txt


#################################### 3.3.3 RNA-RNA Interaction simulation ###########################
python ${scriptdir}/generateSim.py -i $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.genome.interaction.interMolecular.bedpe.gz -o $OUTPUT_DIRECTORY6/sim/$OUTPUT_PREFIX.total.sim -n ${Nsim} > $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.total.sim.log
for pp in $(seq 100); do pigz -p ${THREADS} $OUTPUT_DIRECTORY6/sim/$OUTPUT_PREFIX.total.sim.sim.$pp.tmp.bedpe; done

for pp in $(seq 100)
do
    python ${scriptdir}/getInteraction.py -i $OUTPUT_DIRECTORY6/sim/$OUTPUT_PREFIX.total.sim.sim.${pp}.tmp.bedpe.gz -o $OUTPUT_DIRECTORY6/sim/$OUTPUT_PREFIX.total.sim.sim.${pp}.tmp.bedpe.interaction.txt -l ${genome_info}
    pigz -p ${THREADS} $OUTPUT_DIRECTORY6/sim/$OUTPUT_PREFIX.total.sim.sim.${pp}.tmp.bedpe.interaction.txt
done


#################################### 3.3.4  Indentified Cluster P-value Calculation #################
python ${scriptdir}/calpvalue.py -i $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.total.bedpe.interaction.txt.gz -ip $OUTPUT_DIRECTORY6/sim/$OUTPUT_PREFIX.total.sim.sim -is tmp.bedpe.interaction.txt.gz -n $Nsim -o $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.total.rna-rna.interaction.withpvalue.txt
awk '$21+0.0 <= 0.05' $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.total.rna-rna.interaction.withpvalue.txt > $OUTPUT_DIRECTORY6/$OUTPUT_PREFIX.total.rna-rna.interaction.FDR.txt


