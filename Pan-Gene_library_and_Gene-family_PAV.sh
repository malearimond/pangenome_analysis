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






#Map sequences against reference create paf
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5 /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta unaligned_contigs_renamed.fasta > conts_against_ref.paf "

#remove all sequences from NRR set which coverage <0.6
python paf_drop_coverage.py drop conts_against_ref.paf 0.6 > contigs_to_remove.list
python faSomeRecords.py --fasta unaln_conts.fasta --list contigs_to_remove.list --outfile contigs_dont_match.fasta --exclude

# Cluster all leftover NRR to avoid redundancies with gclust 

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

# Grep sequences from creaeted clusters 
python clu2fa.py contigs_dont_match.sorted.fasta contigs_dont_match.sorted.gclust.clu contigs_dont_match.gclust.fasta

# Check for Contamination: NOTE - maybe next time before clustering?
#!/bin/bash
export BLASTDB="/opt/share/blastdb/ncbi"
blastn -db nt -query ../contigs_dont_match.gclust.masked.fasta -outfmt "6 qseqid bitscore std sscinames stitle staxids sskingdoms" -out blast_hits_unaln_conts_second.txt -num_threads 20 -max_target_seqs 1 -max_hsps 1 -evalue 1e-25
bsub -q bigmem -R "rusage[mem=350000]" -M 450000 "bash blast.sh"

# grep all Aphids found 
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



# FINISHED PAN-GENE LIBRARY
# START CREATING GENE(-FAMILY) PANGENOME 
# First: Annotation ab initio with Helixer, Augustus and SNAP 
# Augustus 
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "augustus --species=arabidopsis contigs_dont_match.gclust.masked.fasta -gff3=on > augustus_novel_genes_masked.gff3"

# Helixer
#!/bin/bash -l
#SBATCH -o ./helixer_%j.out
#SBATCH -e ./helixer_%j.err
#SBATCH -D ./
#SBATCH -J helixer
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=18
#SBATCH --mem=125000
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL,END,TIME_LIMIT_90
#SBATCH --mail-user=jzohren@mpipz.mpg.de
#SBATCH --time=24:00:00

echo 'JOB STARTED'
date


ml purge
ml load intel/21.2.0 impi/2021.2 cuda/11.2
ml singularity
ml

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

OUTDIR=/raven/u/jzohren/230920
cd $OUTDIR
pwd -P
FASTA=$OUTDIR/contigs_dont_match.gclust.fasta

SIF=/raven/u/jzohren/software/helixer-docker_helixer_v0.3.0a0_cuda_11.2.0-cudnn8.sif

echo "running: singularity run --nv --bind $OUTDIR $SIF Helixer.py --fasta-path $FASTA --lineage land_plant --gff-output-path $OUTDIR/Arabis_alpina.MPIPZ.version_5.1_helixer.gff3"
singularity run --nv --bind $OUTDIR $SIF Helixer.py --fasta-path $FASTA --lineage land_plant --gff-output-path $OUTDIR/contigs_dont_match.gclust.gff3


date

# SNAP 


# Grep sequences from GFF with AGAT
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load mambaforge/self-managed/v23.1.0-4
conda activate bases
agat_sp_compare_two_annotations.pl -gff1 augustus_novel_genes_masked.gff3 -gff2 contigs_dont_match.gclust.masked.gff3 -o augustus_helixer_compatison.txt &
agat_sp_compare_two_annotations.pl -gff1 contigs_dont_match.gclust.masked.gff3 -gff2 snap_unaligned_contigs.gff -o snap_helixer_comoparison.txt &
agat_sp_merge_annotations.pl -gff contigs_dont_match.gclust.masked.gff3 -gff snap_unaligned_contigs.gff -gff augustus_novel_genes_masked.gff3 -o merged_annotation.gff

# Merge annotations
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "agat_sp_merge_annotations.pl -gff contigs_dont_match.gclust.masked.gff3 -gff /biodata/dep_coupland/common/Arabis_alpina_resource/annotation/Arabis_alpina_mpipz_v5.1_annotation.gff -out merged_novel_existing_genes_Arabis_alpina.gff"


#####GET TRANSCRIPTS FROM MASKED GCLUST FILE#######
agat_sp_extract_sequences.pl -g contigs_dont_match.gclust.masked.gff3 -f contigs_dont_match.gclust.masked.fasta -o transcripts_helixer.clustered.fasta

# get only complete open reading frames with Transdecoder
TransDecoder.LongOrfs -t transcripts_helixer.clustered.fasta
TransDecoder.Predict -t transcripts_helixer.clustered.fasta


##### MERGE WITH Protein + transcript of orginal Reference
cat Arabis_alpina_mpipz_v5.1_transcripts.fasta >>  pangenome_transcripts_library_final.fasta
cat transcripts_helixer.clustered.fasta >>  pangenome_transcripts_library_final.fasta


##### CLUSTER gene library and grep REPRESENTATIVE SEQUENCES#####
cd-hit -i ../pangenome_protein_library.fna -c 0.98 -o pangenome_protein_library_final.fna
cd-hit-est -i ../pangenome_transcript_library.fasta -c 0.98 -o pangenome_transcripts_library_final.fasta


# Map raw reads against pan-gene library with minimap2

#MAP reads to transcript library:

#!/bin/bash
input_file_SC="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do

echo "Submit ${col1}"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fa ${col2} --secondary=no \ | samtools sort -m 1G -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/bam_files/${col1}_pangenome.bam"
done < "$input_file_SC"
###### OR ########
ID="SAMPLE_ID"
Reference="pangenome_transcripts_library_final.fasta"
raw_reads="RAW_READS"

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "minimap2 -t 20 -ax map-hifi ${reference} ${raw_reads} --secondary=no \ | samtools sort -m 1G -o ${ID}_reads_against_gene_library.bam "

# Calculate Depth 
bedtools genomecov -bga -split -ibam ${ID}_reads_against_gene_library.bam > ${ID}_pan-gene.coverage

# Get all the genes and calculate coverage per gene 
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
while read -r col1 col2 _; do
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/gene_pangenome/
output="/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/gene_pangenome/${col1}_pan-gene.sh"
echo -e "GENE\tSTART\tSTOP\tDEPTH" > ${col1}_pan-gene.txt  && cat ${ID}_pan-gene.coverage  >> ${col1}_pan-gene.txt
echo "
Rscript /netscratch/dep_coupland/grp_fulgione/male/assemblies/Skripte/pangenome/gene_presence_calculator.R.R \"/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/gene_pangenome/${col1}_pan-gene.txt\" " > $output
bsub -q multicore20 -R "rusage[mem=40000]" -M 55000 "bash $output"
done < "$input_file"

# All txt files in same directory and then create (presence_absence_matrix_0.85.tsv) wtih PAV_from_txt.py
python PAV_from_txt.py

#Analyse further in R with Gene-(family) Pangenome Analysis.R 


# NOW GENE FAMILY ANALYSIS WITH ORTHOFINDER
# first get protein sequences depending on txt gene files
#!/bin/bash

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do

rm -r /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences/
mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences/
output="/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences/${col1}_proteins.fasta"

filename="${col1}_present_genes_list.txt"

while IFS= read -r line
do
  echo "${col1}"
  grep -A 1  "$line" /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/pangenome_protein_library.fasta >> $output
done < "$filename"
done < "$input_file"

# RUN orthofinder with A. alpina v5.1, C. rubella, A. lyrate, A. thaliana proteomes (same directory)
https://phytozome-next.jgi.doe.gov/info/Alyrata_v2_1
bsub -q bigmem -R "rusage[mem=450000]" -M 550000 "orthofinder.py -S diamond_ultra_sens -f /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences_0.85/ -t 200"

# Run Assign_orthogroups_to_novel_genes.R
# And use gene_count.tsv for additional R Analysis 
#########################################################



