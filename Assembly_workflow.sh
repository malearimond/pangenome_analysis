#####ASSEMBLY####
#Start: 19.06.2023

#clusters/nodes
marimond@dell-node-12
marimond@hpc001.mpipz.mpg.de

#21 genome assemblies

####################################################################################################################################################
# 1) Heterozygostiy/Genome Size Estimation - Jellyfish
ID=""
path_raw_LR_data=""

jellyfish count -m 21 -s 100M -t 10 -C -F 2 <(zcat ${path_raw_LR_data}) -o ${ID}_mer_counts
#    -m: kmer length, 21 is commonly used
#    -s: size of hash table: should be genome size + extra kmers from seq errors
#    -t: number of threads
jellyfish histo -t 10 --high=1000000 ${ID}_mer_counts >${ID}_reads.histo			#Upload in GenomeScope
jellyfish dump mer_counts.jf > ${ID}_jf_dumps.fa

#with shortreads
sr1=""
sr1=""

jellyfish count -m 21 -s 100M -Q ? -t 10 -C <(zcat ${sr1}) <(zcat ${sr2}) -o ${ID}_sr
jellyfish histo -t 10 --high=1000000 ${ID}_sr > ${ID}_sr.histo 

##################################################################
# 2) de novo assembly using Hifiasm, Flye, Smartdenovo
# 2.1) Hifiasm v0.19
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load hifiasm/v0.19

#without ONT
table_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/assemblies_for_script.txt"
while IFS=$'\t' read -r col1 col2
do

ID="${col1}"
raw_data=${col2}"

#samples with heterozygosity > 0.1:

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/hifiasm_0.19 
bsub -q multicore40 -R rusage[mem=60000] -M 70000 "hifiasm -o ${ID}_hifiasm_0.19.asm -t 40 ${raw_data}"

#samples with heterozygosity < 0.1:

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/hifiasm_0.19 
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ${ID}_hifiasm_0.19.asm -t 40 -l0 ${raw_data}"

#after assembly transform gfa to fasta
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.hap1.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${ID}_hifiasm_0.19.asm.bp.hap2.p_ctg.gfa >  ${ID}_hifiasm_0.19.asm.bp.hap2.p_ctg.fa



done < "$table_file"

#with ONT included

ID=""
ONT_data="PATH_TO_ONT"
raw_data="PATH_TO_RAW_DATA"

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/hifiasm_0.19_ONT/
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ${ID}_hifiasm_0.19_ONT.asm -l0 --ul ${ONT_data} ${raw_data}"

# 2.2) Flye
flye --version        #2.9-b1768 (19.06.2023)

table_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/assemblies_for_script.txt"
while IFS=$'\t' read -r col1 col2
do

ID="${col1}"
raw_data="${col2}"
genome_size="ESTIMATED_GENOME_SIZE"

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/flye/
bsub -q bigmem -R "rusage[mem=200000]" -M 250000 "flye --pacbio-hifi ${raw_data} -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/flye -g ${genome_size} -t 45"
done < "$table_file"


# 2.3) Smartdenovo
table_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/assemblies_for_script.txt"
while IFS=$'\t' read -r col1 col2
do

ID="${col1}"
raw_data="${col2}"

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/smartdenovo/
/netscratch/dep_coupland/grp_fulgione/male/smartdenovo/smartdenovo.pl -t 20 -c 1 -p ${ID} ${raw_data}  > ${ID}_hifi_smdn.mak
bsub -q multicore20 -R "rusage[mem=20000]" -M 35000 "make -f ${ID}_hifi_smdn.mak"

done < "$table_file"

##############################################################################################################################################################################################################
# 3) Purging of haplotigs with purge_dups
#Purge_dups

raw="path/to/data"
ID="ID"
assembly="assembly.fasta"
tool="flye"

#Map reads to assembly and calculate depth +
minimap2 -xasm20 ${assembly} ${raw}  | gzip -c - > ${col1}_${tool}_pb.paf.gz
pbcstat ${ID}_${tool}_pb.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

#split assembly and do self-alignment
split_fa ${assembly} > ${assembly}.split
minimap2 -xasm5 -DP  ${assembly}.split  ${assembly}.split | gzip -c - >  ${assembly}.split.self.paf.gz

#purge haplotigs and overlaps
purge_dups -2 -T cutoffs -c PB.base.cov ${assembly}.split.self.paf.gz > dups.bed 2> purge_dups.log

#get purged sequences and haplotigs
get_seqs -dups.bed ${assembly}


bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bash purge_dups.sh"


#################################################################################################################################
# 4) Phasing of Assemblies with HapDup 
ID=""
tool="ASSEMBLY_TOOL"
raw_data="path_to_raw_data"
assembly="path_to_assembly"

#create bam file between assembly and raw reads with minimap2
minimap2 -t 20 -ax map-hifi ${assembly} ${raw_data} --secondary=no \ | samtools sort -m 1G -o ${ID}_${tool}_pb_aligned.bam
samtools index -b ${ID}_${tool}_pb_aligned.bam
#run HapDup on bam file
singularity pull docker://mkolmogo/hapdup:0.12
HD_DIR=`pwd`

singularity exec --bind $HD_DIR hapdup_0.12.sif \
  hapdup --assembly ${assembly} \
   --bam ${ID}_${tool}_pb_aligned.bam \
   --bam-index ${ID}_${tool}_pb_aligned.bam.bai \
   --out-dir ${ID}_hapdup_out -t 10 --rtype hifi

##################################################################################################################################
# 5) Merge assembly sets with quickmerge
#nucmer 4.0.0beta2 # Nucmer aligns the two assemblies so that the merger can find the correct splice sites:
# assembly with the longest contiguity (N50) selected as query, the assembly with second longest N50 was used as reference to join query assembly
#The resulting assembly was further improved by using the third assembly as reference


#####run  EXPORT#########
query="query_assembly"
reference="reference_assembly"

export PATH=/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge/MUMmer3.23:$PATH
merge_wrapper.py ${reference} ${query}

############################################################################################
# 5) Scaffolding with RagTag
reference="/netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta"
assembly="ASSEMBLY"
ID=""
mapper="PATH_TO_MINIMAP_OR_NUCMER"

export PATH="/home/marimond/.conda/envs/bases/bin/:$PATH"
# without chimeric contig correction
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold ${reference} ${assembly} -o ${ID}_RagTag -C --aligner ${mapper} -t 20"
# with chimeric contig correction
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py correct ${reference} ${assembly}  -T corr -o ${ID}_RagTag --aligner ${mapper} -t 20"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold ${reference} ${assembly} -o ${ID}_RagTag -C --aligner ${mapper} -t 20"


export PATH="/home/marimond/.conda/envs/bases/bin/:$PATH"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/SI01_005/merged_flye_hifiasm/merged_out.fasta -C -t 20"

####################################################
# 6) Gapclosing with LR Gapcloser
export PATH="/netscratch/dep_coupland/grp_fulgione/male/assemblies/LR_Gapcloser/src/:$PATH"
bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "bash LR_Gapcloser.sh -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag_nucmer/ragtag_wo_corr/ragtag.scaffold.fasta -l SE02_002_raw.fasta -o LR_Gapcloser/"

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
lsfenv
while read -r col1 col2 _; do
  samtools faidx /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/ragtag.scaffold.fasta
  echo "Index of ${col1} created"
  seqtk seq -a ${col2} > /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/${col1}_raw.fasta
  export PATH="/netscratch/dep_coupland/grp_fulgione/male/assemblies/LR_Gapcloser/src/:$PATH"
  bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "bash LR_Gapcloser.sh -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/ragtag.scaffold.fasta -l /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/${col1}_raw.fasta -m 1000 -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser"
done < "$input_file"

echo "all gapclosing jobs sended"

#################################################################################
# 7) Evaluation of assemblies
# 7.1) Create BAM and PAF

#CREATE PAF AND BAM for Evaluation

#!/bin/bash
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do
  echo "Submit ${col1}"
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly/
  bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi AT03_004_final_assembly_minimap.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078_210408_170302.Q20.fastq.gz --secondary=no \ | samtools sort -m 1G -o AT03_004_reads_against_assembly.bam"
  bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5 /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta ${col1}_assembly_only_chr.fasta > ${col1}_scaffolds_against_ref.paf"

done < "$input_file"

#create coverage plots 
while read -r col1 col2 _; do
  output_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly/long_read_depth_plots.sh"
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly/
  echo "#!/bin/bash" >> $output_file
  echo "
  samtools depth -a ${col1}_reads_against_assembly.bam > ${col1}_LR.depth
  awk '{print $1}' ${col1}_LR.depth | uniq
  awk '$1 == \"chr1_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr1.coverage
  awk '$1 == \"chr2_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr2.coverage
  awk '$1 == \"chr3_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr3.coverage
  awk '$1 == \"chr4_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr4.coverage
  awk '$1 == \"chr5_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr5.coverage
  awk '$1 == \"chr6_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr6.coverage
  awk '$1 == \"chr7_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr7.coverage
  awk '$1 == \"chr8_RagTag\" {print $0}' ${col1}_LR.depth > LR_chr8.coverage
  awk '$1 == \"Chloroplast_Aa_NC_023367.1_RagTag\" {print $0}' ${col1}_LR.depth > LR_cloroplast.coverage
  awk '$1 == \"Mitochondrium_Aa_NC_037070.1_RagTag\" {print $0}' ${col1}_LR.depth > LR_mitochondrium.coverage

  Rscript ${col1}_temp_R_script.R \"${col1}\"" >> "$output_file"
  bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bash $output_file"

done < "$input_file"

#For PAF files run PAF_Skript.R


# 7.2) Run BUSCO and QUAST
# 3) Evaluation of assemblies
# 3.1) QUAST
ID=""

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/flye
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/flye/assembly.fasta -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/flye/quast -e -f -b 
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/smartdenovo
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/smartdenovo/${ID}_hifi.dmo.cns -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/smartdenovo/quast -e -f -b 

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/hifiasm
/home/marimond/.conda/envs/bases/bin/quast.py ${ID}_hifiasm_0.19.asm.bp.hap1.p_ctg.fa -o hap1_quast -e -f -b 
/home/marimond/.conda/envs/bases/bin/quast.py ${ID}_hifiasm_0.19.asm.bp.hap2.p_ctg.fa -o hap2_quast -e -f -b 
/home/marimond/.conda/envs/bases/bin/quast.py ${ID}_hifiasm_0.19.asm.bp.p_ctg.fa -o both_quast -e -f -b 

# 3.2) BUSCO 
assembly="ASSEMBLY"
busco -i ${assembly} -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_${assembly} -m genome -f 
busco -i SE02_002_hifiasm_0.19.asm.bp.hap2.p_ctg.fa -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_hap2 -m genome -f &
busco -i ES20_028_final_assembly_minimap.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_ES24_001 -m genome -e

########### FINAL EVALUATION ##################
#!/bin/bash

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/Homozygots_only_LR.txt"
analysis="/netscratch/dep_coupland/grp_fulgione/male/assemblies/analysis_final_assemblies_quast.txt"
rm $analysis
while read -r col1 col2 _; do

echo "BUSCO" >> $analysis
 grep -B 1 'BUSCOs' /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly_minimap/busco_${col1}/run_brassicales_odb10/short_summary.txt >> $analysis


echo "Quast" >> $analysis
  cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly_minimap/quast_${col1}/report.txt >> $analysis
  target_character="N"
  fasta_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly_minimap/${col1}_final_assembly_minimap.fasta"
  grep -v "^>" "$fasta_file" | tr -d '\n' > sequences.txt
  awk -v char="$target_character" '{print gsub(char, "")}' sequences.txt | paste -sd+ | bc >> $analysis
  rm sequences.txt
  samtools faidx $fasta_file
  echo "chromosome lenghts:" >> $analysis
  cut -f1,2 ${fasta_file}.fai >> $analysis
  echo "total length:" >> $analysis
  awk -F' ' '{sum+=$2} END {print sum}' ${fasta_file}.fai >> $analysis
  cut -f2,3,4 /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_minimap/ragtag_output/ragtag.scaffold.confidence.txt | awk '{ total_grouping += $1; total_location += $2; total_orientation += $3; count++ } END { print "Mean Grouping Confidence:", total_grouping/count; print "Mean Location Confidence:", total_location/count; print "Mean Orientation Confidence:", total_orientation/count }' >> $analysis
  echo "--------------------------------------------------------------------------------------------" >> $analysis
done < "$input_file"

###############################################################
#START PAN GENE LIBRARY 


cd-hit -i transcripts_helixer.clustered.fasta.transdecoder.pep -c 0.98 -o transcripts_helixer.clustered.fasta.transdecoder.fna

bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi SI01_005_assembly_only_chr.fasta /netscratch/dep_coupland/grp_fulgione/danijel/globus/891_1_CELL1/m54274Ue_230705_214109.bc1001--bc1001.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o SI01_005_reads_against_assembly.bam"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5 /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta SI01_005_assembly_only_chr.fasta > SI01_005_scaffolds_against_ref.paf"


#####Gene Liftover#####
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load python/v3.8.0

/opt/share/software/scs/appStore/selfStorage/marimond/stretch/mambaforge/v23.1.0-4/envs/liftoff/bin/liftoff
marimond@dell-node-12:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/liftoff$ /opt/share/software/scs/appStore/selfStorage/marimond/stretch/mambaforge/v23.1.0-4/envs/liftoff/bin/liftoff -g ../../Reference_v5.1_analysis/Arabis_alpina_mpipz_v5.1_annotation.gff ../final_assembly/SE02_002_final_assembly.fasta ../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
lsfenv
while read -r col1 col2 _; do
  mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' ../pilon/pilon.fasta | tr "\t" "\n" > /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly/${col1}_final_assembly.fasta
  grep -A 1 'chr' ${col1}_final_assembly.fasta > ${col1}_assembly_only_chr.fasta
  bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "/opt/share/software/scs/appStore/selfStorage/marimond/stretch/mambaforge/v23.1.0-4/envs/liftoff/bin/liftoff -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_only_chr.gff -o ${col1}_annotation.gff -flank 0.2 -copies ${col1}_assembly_only_chr.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated_only_chr.fasta"
done < "$input_file"

bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "/opt/share/software/scs/appStore/selfStorage/marimond/stretch/mambaforge/v23.1.0-4/envs/liftoff/bin/liftoff -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_only_chr.gff -o NO01_005_annotation.gff -flank 0.2 -copies NO01_005_assembly_only_chr.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated_only_chr.fasta"

#########BUSCI QUAST final assembly################
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/Heterozygots.txt"
while read -r col1 col2 _; do
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly_minimap_purged
  samtools faidx ${col1}_final_assembly_minimap.fasta
  bsub -q multicore20 -R "rusage[mem=20000]" -M 35000 "/home/marimond/.conda/envs/bases/bin/quast.py ${col1}_final_assembly_minimap.fasta  -o quast_purged_ass_${col1} -e -f -b"
  bsub -q multicore20 -R "rusage[mem=20000]" -M 35000 "busco -i ${col1}_final_assembly_minimap.fasta  -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_purged_ass_${col1} -m genome -f"
done < "$input_file"
busco -i ES06_006_final_assembly.fasta  -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_merged_ES06_006 -m genome -f


source /opt/share/software/scs/appStore/modules/init/profile.sh
module load mambaforge/self-managed/v23.1.0-4





###################### PanGenome ######################
#!/bin/bash
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
#improve transcript databank of unknown sequences
#look at all unaligned contigs of Chr0
mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/
grep 'Chr0' /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/ragtag.scaffold.agp | grep 'contig_' | awk '{print $6}' > ${col1}_unaln.contig.txt

#Extract
# Define input and output files
contig_list_file="${col1}_unaln.contig.txt"
fasta_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/merged_flye_hifiasm/merged_out.fasta"
output_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/unaligned_contigs.fasta"
final_output="/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/unaligned_contigs_renamed.fasta"

# Create an associative array to store contigs as keys
declare -A contigs

# Read the list of contigs into the associative array
while IFS= read -r contig; do
    contigs["$contig"]=1
done < "$contig_list_file"

# Initialize a flag variable to track whether a contig is found
found=0



while read -r col1 col2 _; do
sed -e 's/contig/${col1}_contig/g' $output_file >> $final_output
echo "${col1} contigs extracted"
head $final_output
#merged
done < "$input_file"
#check for redundancies
#annotate CD-Hit
#merge gff




export PATH=/netscratch/dep_coupland/grp_fulgione/male/assemblies/EUPAN-v0.44/bin:$PATH

export LD_LIBRARY_PATH=/netscratch/dep_coupland/grp_fulgione/male/assemblies/EUPAN-v0.44/lib/:$LD_LIBRARY_PATH
export PERL5LIB=/netscratch/dep_coupland/grp_fulgione/male/assemblies/EUPAN-v0.44/lib/:$PERL5LIB
source /netscratch/dep_coupland/grp_fulgione/male/assemblies/EUPAN-v0.44/bin/eupan_cmd.sh



export PATH=$PATH:/netscratch/dep_coupland/grp_fulgione/male/assemblies/HUPAN-v1.04/bin:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/netscratch/dep_coupland/grp_fulgione/male/assemblies/HUPAN-v1.04/lib/:
export PERL5LIB=$PERL5LIB:/netscratch/dep_coupland/grp_fulgione/male/assemblies/HUPAN-v1.04/lib/:
source /netscratch/dep_coupland/grp_fulgione/male/assemblies/HUPAN-v1.04/hupan_cmd.sh


hupanSLURM alignContig -s merged_out.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/data_contigs /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaligned_contigs /opt/share/software/bin/mummer /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta

#Split sets
import random

def split_fasta(input_file, output_prefix, num_datasets):
    # Read the input FASTA file into a dictionary
    sequences = {}
    current_sequence = ''
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_sequence = line.strip()
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line.strip()

    # Shuffle the sequences randomly
    sequence_keys = list(sequences.keys())
    random.shuffle(sequence_keys)

    # Split the sequences into smaller datasets
    chunk_size = len(sequence_keys) // num_datasets
    for i in range(num_datasets):
        dataset_name = '{}_{}.fasta'.format(output_prefix, i)
        with open(dataset_name, 'w') as f:
            for key in sequence_keys[i * chunk_size:(i + 1) * chunk_size]:
                f.write(key + '\n')
                f.write(sequences[key] + '\n')

if __name__ == '__main__':
    input_file = 'input.fasta'  # Replace with your input FASTA file
    output_prefix = 'dataset'  # Prefix for output dataset filenames
    num_datasets = 5  # Number of smaller datasets to create

    split_fasta(input_file, output_prefix, num_datasets)


##CDHIT remove redundancies
#Map sequences against reference

#!/bin/bash


input="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
output="/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/unaln_conts.fasta"


while read -r col1 col2 _; do
grep -A 1 -f /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/contigs_per_sample/${col1}_unaln.contig.txt /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/merged_flye_hifiasm/merged_out.fasta | grep -v '^--$' | sed s/contig/contig_${col1}/ >> $output

done < $input


bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax asm5 ../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta unaln_conts.fasta  --secondary=no \ | samtools sort -m 1G -o unaln_conts_sort.sam"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "blastn -db /opt/share/blastdb/ncbi/nt -query contigs_dont_match.gclust.fasta -outfmt 6 -out blast_hits_nt -num-threads 20"

python3 getTax.py blastnt/unalign_clsutered.blastnt 2 ~/data/db/taxonomy_20200715/nucl_gb.accession2taxid ~/data/db/taxonomy_20200715/rankedlineage.dmp > blastnt/unalign_clsutered.tax
python3 filter_exclude_contig.py blastnt/unalign_clsutered.blastnt blastnt/uunalign_clsutered.tax 0 0 0 > blastnt/unalign_clsutered.polution
~/tool/faSomeRecords ${query} -exclude blastnt/unalign_clsutered.polution blastnt/unalign_clsutered_clean.fa
https://github.com/SJTU-CGM/TGSRICEPAN/blob/master/Genomes2Pan.sh

#!/bin/bash
export BLASTDB="/opt/share/blastdb/ncbi"
blastn -db nt -query ../contigs_dont_match.gclust.masked.fasta -outfmt "6 qseqid bitscore std sscinames stitle staxids sskingdoms" -out blast_hits_unaln_conts_second.txt -num_threads 20 -max_target_seqs 1 -max_hsps 1 -evalue 1e-25
bsub -q bigmem -R "rusage[mem=350000]" -M 450000 "bash blast.sh"


#!/bin/bash
grep 'Myzus' blast_hits_unaln_conts_second.txt | less -S > Myzus_contaminations.txt
grep 'Bacteria' blast_hits_unaln_conts_second.txt | less -S >> Myzus_contaminations.txt
grep 'Periphyllus' blast_hits_unaln_conts_second.txt | less -S > Periphyllus_contamination.txt
grep 'gossy' blast_hits_unaln_conts_second.txt | less -S > Aphis_gossypii_contamination.txt
grep 'Acyrthosiphon' blast_hits_unaln_conts_second.txt | less -S > Acyrthosiphon_contaminationo.txt
grep 'Melanaphis sacchari' blast_hits_unaln_conts_second.txt | less -S > Melanaphis_sacchari_contaminationo.txt
rm all_contigs_contaminated.txt
awk '{print $1}' Acyrthosiphon_contaminationo.txt >> all_contigs_contaminated.txt
awk '{print $1}' Aphis_gossypii_contamination.txt >> all_contigs_contaminated.txt
awk '{print $1}' Myzus_contaminations.txt >> all_contigs_contaminated.txt
awk '{print $1}' Periphyllus_contamination.txt >> all_contigs_contaminated.txt
awk '{print $1}' Melanaphis_sacchari_contaminationo.txt >> all_contigs_contaminated.txt
sort all_contigs_contaminated.txt | uniq > all_contigs_contaminated_unique.txt
echo "in total:"
wc -l all_contigs_contaminated_unique.txt
echo "contaminated conigs"

#REMOVE Myzus_contaminations
bash /netscratch/dep_coupland/grp_fulgione/male/assemblies/Skripte/pangenome/remove_contaminations.sh
##ANALYSIS
# set path to local BLAST databases
export BLASTDB="/opt/share/blastdb/ncbi"
blobtools create -i ../contigs_dont_match.gclust.masked.fasta -t blast_hits_unaln_conts.txt --nodes /opt/share/blastdb/taxonomy/nodes.dmp --names /opt/share/blastdb/taxonomy/names.dmp -o unaligned_contigs_TaxIDs

bsub -q bigmem -R "rusage[mem=400000]" -M 430000 "bash /netscratch/dep_coupland/grp_fulgione/male/assemblies/Skripte/analysis_final_assemblies/blobtools.sh"

#Filter out unmapped sequences
samtools view -F 4 -b ref_cont_sort.sam > mapped_reads.bam
samtools bam2fq -c 0 unmapped_reads.bam > unmapped_reads.fq




python paf_drop_coverage.py drop conts_against_ref.paf 0.6 > contigs_to_remove.list

python faSomeRecords.py --fasta unaln_conts.fasta --list contigs_to_remove.list --outfile contigs_dont_match.fasta --exclude

export PATH=/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/gclust/:$PATH
/gclust/script/sortgenome.pl --genomes-file contigs_dont_match.fasta \
	--sortedgenomes-file contigs_dont_match.sorted.fasta

gclust -loadall \
	-minlen 20 \
	-both \
	-nuc \
	-threads 40 \
	-ext 1 \
	-sparse 2 \
	-memiden 90 \
	contigs_dont_match.sorted.fasta > contigs_dont_match.sorted.gclust.clu

python clu2fa.py contigs_dont_match.sorted.fasta contigs_dont_match.sorted.gclust.clu contigs_dont_match.gclust.fasta
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/unaln_contigs/contigs_dont_match.gclust.fasta
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5 /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta ES18_040_final_assembly.fasta > ES18_040_all_chr_against_ref_ragtag_nucmer.paf"

https://www.plabipd.de/helixer_main.html

bsub -q bigmem -R "rusage[mem=200000]" -M 220000 "blastclust -i contigs_dont_match.gclust.fasta -a 20 -p F -n"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "augustus --species=arabidopsis contigs_dont_match.gclust.fasta -gff3=on > novel_genes.gff3"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "augustus --species=arabidopsis contigs_dont_match.gclust.masked.fasta --gff3=on > augustus_novel_genes_masked.gff3"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "augustus --species=arabidopsis contigs_dont_match.gclust.masked.fasta -gff3=on > augustus_novel_genes_masked.gff3"
bsub -q bigmem -R "rusage[mem=305000]" -M 315000 "blastp -db /opt/share/blastdb/ncbi/nt -query contigs_dont_match.gclust.fasta -outfmt 6 -out blast_hits_nr -num_threads 20 -evalue 1e-5 -best_hit_overhang 0.25 -perc_identity 0.95"
BLASTDB = /opt/share/blastdb/ncbi/FASTA/nt
-outfmt "6 qseqid sseqid slen qstart qend length mismatch gapopen gaps sseq sskingdoms"
export BLASTDB
/biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078_210408_170302.Q20.fastq.gz

n9T9dJ7ppn


source /opt/share/software/scs/appStore/modules/init/profile.sh
module load mambaforge/self-managed/v23.1.0-4
conda activate bases
agat_sp_extract_sequences.pl

agat_sp_compare_two_annotations.pl -gff1 augustus_novel_genes_masked.gff3 -gff2 contigs_dont_match.gclust.masked.gff3 -o augustus_helixer_compatison.txt &
agat_sp_compare_two_annotations.pl -gff1 contigs_dont_match.gclust.masked.gff3 -gff2 snap_unaligned_contigs.gff -o snap_helixer_comoparison.txt &
agat_sp_merge_annotations.pl -gff contigs_dont_match.gclust.masked.gff3 -gff snap_unaligned_contigs.gff -gff augustus_novel_genes_masked.gff3 -o merged_annotation.gff

bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "agat_sp_merge_annotations.pl -gff contigs_dont_match.gclust.masked.gff3 -gff /biodata/dep_coupland/common/Arabis_alpina_resource/annotation/Arabis_alpina_mpipz_v5.1_annotation.gff -out merged_novel_existing_genes_Arabis_alpina.gff"


#####GET TRANSCRIPTS FROM MASKED GCLUST FILE#######
agat_sp_extract_sequences.pl -g contigs_dont_match.gclust.masked.gff3 -f contigs_dont_match.gclust.masked.fasta -o transcripts_helixer.clustered.fasta

#####MERGE TRASNCRIPTS WITH REFERENCE######
TransDecoder.LongOrfs -t transcripts_helixer.clustered.fasta
TransDecoder.Predict -t transcripts_helixer.clustered.fasta
source /netscratch/common/Saurabh/envMod/init/bash
module load interproscan/v5.52-86.0

interproscan.sh -i pangenome_final_proteins.fna --excl-applications ProSiteProfiles,ProSitePatterns --goterms --formats TSV --output-file-base out --cpu 4

#takes a proteome fasta file and splits into 100 parts
#run on terminal without lsf
#load the seqkit module
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load seqkit/v0.14.0
#use seqkit to split the fasta file in 100 parts and keep the output in parts/ folder
seqkit split -p 100 pangenome_final_proteins.fna -O parts/
#this will generate 100 files in the parts/ folder named  carHir1_lprot.part_001.fa  carHir1_lprot.part_002.fa
#sdf

##### MERGE WITH Protein + transcript of orginal Reference#####


#####CLUSTER: TAKE REPRESENTATIVE SEQUENCES#####
cd-hit -i ../pangenome_protein_library.fna -c 0.98 -o pangenome_protein_library_final.fna

cd-hit-est -i ../pangenome_transcript_library.fasta -c 0.98 -o pangenome_transcripts_library_final.fasta

###### OrthoFinder to get homologous seqeunces/gene families####

bsub -q bigmem -R "rusage[mem=450000]" -M 550000 "orthofinder.py -S diamond_ultra_sens -f /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences_0.85/ -t 200"
#Photozyme Alyrata v2.1 381   https://phytozome-next.jgi.doe.gov/info/Alyrata_v2_1
Exonerate
v2.2.0102 399 (--percent 70 --minintron 10 --maxintron 60000)
#####CLEAN ORTHOFINDER FILES##########
awk -F'\t' '$2 ~ /^_contig/' Orthogroups.tsv > Orthogroups_novel_genes_old.tsv
sed -i -E 's/ORF type_complete len_|_+__score=|score=[0-9]+\.?[0-9]* //' your_file.txt
sed -i 's/_contig_[^[:space:]]*//g' Orthogroups_novel_genes.tsv
sed -i 's/ORF type_complete len_[0-9]+\.?[0-9]*//g; s/_+__score=//g' Orthogroups_novel_genes.tsv
sed -i 's/_+__score=[0-9][0-9]//g' Orthogroups_novel_genes.tsv
sed -i 's/len_[0-9][0-9][0-9]//g' Orthogroups_novel_genes.tsv
sed -i 's/GENE.  ORF type_complete//g' Orthogroups_novel_genes.tsv
sed -i 's/GENE.  ORF type_5prime_partial//g' Orthogroups_novel_genes.tsv
sed -i 's/len_//g' Orthogroups_novel_genes.tsv
sed -i 's/ORF type_3prime_partial//g' Orthogroups_novel_genes.tsv
sed -i 's/GENE.//g' Orthogroups_novel_genes.tsv
sed -i 's/[0-9].[0-9][0-9]*//g' Orthogroups_novel_genes.tsv


#######PANGENOME#######

minimap2 -a reference_transcripts.fasta raw_long_reads.fastq > alignments.sam
samtools view -b -o alignments.bam alignments.sam
samtools sort -o alignments_sorted.bam alignments.bam
samtools index alignments_sorted.bam
genomecov -bga -split


/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/EUPAN-v0.44/bin/


# Calculate the threshold (e.g., 95% of total length)
total_length=1000000  # Replace with your actual total length
threshold=$(echo "$total_length * 0.95" | bc)
your_coverage_file="DE01_001_pan-gene.coverage"
# Filter hits with coverage greater than or equal to the threshold
awk -v threshold="$threshold" '($4 >= threshold)' your_coverage_file | cut -f1 > high_coverage_hits.txt


input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
output="/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences/counts_proteins/protein_per_sample.txt"
while read -r col1 col2 _; do
  echo "${col1} has:" >> $output
  grep '>' /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences/${col1}_proteins.fasta | wc -l >> $output
  echo "protein CDS" >> $output
  echo "------------------------------" >> $output
done < "$input_file"


bedtools genomecov -bga -split -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fasta -ibam ES24_001_splitted_pangenome.bam > ES24_001_splitted.coverage
awk '($4 > 3)' CZ01_005_pan-gene.coverage | cut -f1 | sort | uniq > SL01_005_genes_present_2.txt
bedtools genomecov -bga -split -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fasta -ibam AT03_004_splitted_pangenome.bam > AT04_splitted.coverage
bedtools genomecov -bga -split -ibam CZ01_005_reads_against_reference.bam > CZ01_005_reads_against_reference.coverage
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES27_024/
bedtools genomecov -bga -split -ibam ES27_024_reads_against_reference.bam > ES27_024_reads_against_reference.coverage
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/
bedtools genomecov -bga -split -ibam AT03_004_reads_against_reference.bam > AT03_004_reads_against_reference.coverage
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES04_014/
bedtools genomecov -bga -split -ibam ES04_014_reads_against_reference.bam > ES04_014_reads_against_reference.coverage
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE07_013/
bedtools genomecov -bga -split -ibam SE07_013_reads_against_reference.bam > SE07_013_reads_against_reference.coverage

eupan geneCov /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/ /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/eupan_geneCov/ /biodata/dep_coupland/common/Arabis_alpina_resource/annotation/Arabis_alpina_mpipz_v5.1_annotation.gtf
Arabidopsis thaliana TAIR10
Capsella rubella v1.1



#### Clean contamination #####
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -x map-hifi ES20_028_conts_and_sequences_to_rm_contigs.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210729_005338.hifi_reads.fastq.gz > ES20_028_contamination_free.paf"
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -x map-hifi ES18_conts_and_sequences_to_rm_contigs.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210729_191323.hifi_reads.fastq.gz > ES18_040_contamination_free.paf"

makeblastdb -in ES20_028_conts_and_sequences_to_rm_contigs.fasta -dbtype nucl -out ES20_028_con_free
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "blastn -db ES20_028_con_free -query /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210729_005338.hifi_reads.fastq.gz -outfmt 6 -out ES20_contamination_blast_hits -num_threads 20"

bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "blastn -db ES18_040_con_free -query /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210729_191323.hifi_reads.fastq.gz -outfmt 6 -out ES18_contamination_blast_hits -num_threads 20"


bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_protein_library_final.fa /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210722_072017.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o CZ01_005_pangenome.bam"
output="/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences/CZ01_005_proteins.fasta"

filename="CZ01_005_genes_present_2.txt"

while IFS= read -r line
do
  echo "CZ01_005"
  grep -A 1  "$line" /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/pangenome_protein_library.fasta >> $output
done < "$filename"

############ Gene analysis ##############
filename="genes_from_0_to_1.txt"
output_file="gene_analysis.txt"
while IFS= read -r line
do
grep '$line' merged_novel_existing_genes_Arabis_alpina.gff | grep 'mRNA' >> $output_file
echo $line
done < "$filename"
################# Cutoff differences #######################

filename="gene_cutoff_difference.txt"
output_file="gene_analysis_cutoff_differences.gff"
#output_file_2="gene_analysis_cutoff_differences.fasta"
while IFS= read -r line
do
grep '$line' annotation_files/merged_novel_existing_genes_Arabis_alpina.gff | grep 'mRNA' >> $output_file
echo $line
done < "$filename"



#########################################################
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
while read -r col1 col2 _; do
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/gene_pangenome/
output="/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/gene_pangenome/${col1}_pan-gene.sh"
echo -e "GENE\tSTART\tSTOP\tDEPTH" > NO04_subsampled_pan-gene_splitted.txt  && cat NO04_subsampled_pangenome.coverage  >> NO04_subsampled_pan-gene_splitted.txt
echo "
Rscript /netscratch/dep_coupland/grp_fulgione/male/assemblies/Skripte/pangenome/gene_presence_calculator.R.R \"/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/gene_pangenome/${col1}_pan-gene.txt\" " > $output
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "bash $output"
done < "$input_file"
################################################################################
AB000835.1
esearch -db nucleotide -query "AB000835.1" | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS=", ";}{print $3;}'

Rscript /netscratch/dep_coupland/grp_fulgione/male/assemblies/Skripte/pangenome/gene_presence_calculator.R.R "/netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/gene_pangenome/AT03_004_pan-gene.txt"
Rscript /netscratch/dep_coupland/grp_fulgione/male/assemblies/Skripte/pangenome/gene_presence_calculator.R.R \"/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/gene_pangenome/${col1}_pan-gene.txt\" > $output

/netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/gene_pangenome/

######################## HIST 0 ANALYSIS for SPLITTED ###########
mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/IT19_017/
mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/IT19_017/gene_pangenome/
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/IT19_017/gene_pangenome/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_03_16/demultiplex.bc1010_BAK8A_OA--bc1010_BAK8A_OA.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o IT19_017_gene_pangenome.bam"

mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/IT20_007/
mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/IT20_007/gene_pangenome/
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/IT20_007/gene_pangenome/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_03_16/demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o IT20_007_gene_pangenome.bam"





##########################################################
#COVERAGE WITH REFERENCE
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078_210408_170302.Q20.fastq.gz --secondary=no \ | samtools sort -m 1G -o AT03_004_reads_against_reference.bam"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES04_014/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210726_124221.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o ES04_014_reads_against_reference.bam"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE07_013/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210720_085032.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o SE07_013_reads_against_reference.bam"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/CZ01_005/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210722_072017.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o CZ01_005_reads_against_reference.bam"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES27_024/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210723_131952.hifi_reads.fastq.gz --secondary=no \ | samtools sort -m 1G -o ES27_024_reads_against_reference.bam"

########################### Calculate median depth #################
#!/bin/bash

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
while read -r col1 col2 _; do

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly_minimap/

R_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly_minimap/${col1}_median_depth.R"
echo "
#!/usr/bin/env Rscript

#Load Packages

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
DE01_001_depth <- read.table(name, header =T )
median <- median(DE01_001_depth[,1])
write.table(median, print0("median_", name)" > $R_file

awk '{print $3}' ${col1}_LR.depth > depth_${col1}
bsub -q multicore20 -R "rusage[mem=20000]" -M 35000 "Rscript $R_file depth_${col1}"
done < "$input_file"


#########################################################
#MAP reads to transcript library:

#!/bin/bash
input_file_SC="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do

echo "Submit ${col1}"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fa ${col2} --secondary=no \ | samtools sort -m 1G -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/${col1}_pangenome.bam"
done < "$input_file_SC"

#GET COVERAGE FROM BAM FILES

#!/bin/bash
input_file_SC="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do

mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bedtools genomecov -bga -split -ibam /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/${col1}_pangenome.bam > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/${col1}_pangenome.coverage"

done < "$input_file_SC"

#!/bin/bash
input_file_SC="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do
echo -e "GENE\tSTART\tSTOP\tDEPTH" > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/${col1}_pan-gene.txt  && cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/${col1}_pangenome.coverage >> /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/${col1}_pan-gene.txt
done < "$input_file_SC"

input_file_SC="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do
echo "${col1_}_tmpfile"
grep -v '^--' -i ${col1}_proteins.fasta > ${col1_}_tmpfile && mv ${col1_}_tmpfile ${col1}_proteins.fasta
done < "$input_file_SC"

##########Subsample data############
input_file_SC="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
output="comparison_genes.txt"
while read -r col1 col2 _; do
echo "${col1}" >> $output
grep -v 'contig' ${col1}_present_genes_list_0.85.txt | wc -l >> $output
done < "$input_file_SC"



bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fa NO01_005_split_data_0.fasta --secondary=no \ | samtools sort -m 1G -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/NO01_005_split_data_pangenome.bam"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bedtools genomecov -bga -split -ibam /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/NO01_005_split_data_pangenome.bam > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data.coverage"
echo -e "GENE\tSTART\tSTOP\tDEPTH" > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data_pan-gene.txt  && cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data.coverage >> /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data_pan-gene.txt

bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bedtools genomecov -bga -split -ibam /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/NO04_subsampled_pangenome.bam > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO04_subsampled_pangenome.coverage"
rm /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data_pan-gene.txt
echo -e "GENE\tSTART\tSTOP\tDEPTH" > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data_pan-gene.txt  && cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data.coverage >> /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data_pan-gene.txt



bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fa NO01_005_split_data_0.fasta --secondary=no \ | samtools sort -m 1G -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/NO01_005_split_data_pangenome.bam"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bedtools genomecov -bga -split -ibam /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/NO01_005_split_data_pangenome.bam > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data.coverage"
echo -e "GENE\tSTART\tSTOP\tDEPTH" > /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data_pan-gene.txt  && cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data.coverage >> /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/coverage_files/NO01_005_split_data_pan-gene.txt


##################################



import pysam

import pandas as pd

# List of BAM file paths (replace with your actual file paths)
bam_files = ["assembly1.bam", "assembly2.bam", ..., "assembly20.bam"]

# Initialize a dictionary to store presence/absence information
presence_absence_dict = {}

# Loop through BAM files
for bam_file in bam_files:
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Initialize a set to store the transcripts present in this assembly
        present_transcripts = set()

        # Iterate through alignments
        for alignment in bam:
            # Extract the reference sequence ID (transcript ID)
            transcript_id = alignment.reference_name

            # Update the set to indicate presence
            present_transcripts.add(transcript_id)

        # Store the presence/absence information for this assembly
        presence_absence_dict[bam_file] = present_transcripts

# Create a DataFrame from the presence/absence dictionary
matrix = pd.DataFrame(presence_absence_dict)

# Save the matrix to a file or perform further analysis
matrix.to_csv("presence_absence_matrix.csv", index=False)
