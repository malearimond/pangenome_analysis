#####ASSEMBLY####
#Start: 19.06.2023
marimond@dell-node-12
#path system
/netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/hifiasm/${ID}_hifiasm.p_ctg.fa
/netscratch/dep_coupland/grp_fulgione/male/assemblies/<file ID>/<data_type (pacbio_hifi; ONT)>/<assembler>/<assembly>
${ID}=it20_007
${raw_data}

sed 's/\r//' #to remove ^M
###Data:
#${ID} PacBio ${path_raw_LR_data}

####################################################################################################################################################
###Genome Size Estimation - Jellyfish###
jellyfish count -m 21 -s 100M -t 10 -C -F 2 <(zcat /netscratch/dep_coupland/grp_fulgione/danijel/globus/891_6_CELL_merged/GR01_001_merged.fq.gz)
#    -m: kmer length, 21 is commonly used
#    -s: size of hash table: should be genome size + extra kmers from seq errors
#    -t: number of threads
jellyfish histo -t 10 --high=1000000 mer_counts.jf > GR01_001_reads.histo
jellyfish dump mer_counts.jf > ${col1}_jf_dumps.fa

#shortreads
#nor2
jellyfish count -m 21 -s 100M -Q ? -t 10 -C <(zcat /biodata/dep_coupland/grp_fulgione/rawdata/alpina_24_03_2021/F20FTSEUHT1304_ARAsovR/Clean/233/FP100002218TR_L01_SP2101080192-541_1.fq.gz) <(zcat /biodata/dep_coupland/grp_fulgione/rawdata/alpina_24_03_2021/F20FTSEUHT1304_ARAsovR/Clean/233/FP100002218TR_L01_SP2101080192-541_2.fq.gz) -o nor02_001_sr_qual30 &
jellyfish histo -t 10 --high=1000000 nor02_001_sr_qual30 > nor02_001_sr_qual30.histo &
jellyfish histo -t 10 --high=1000000 nor02_001_sr_qual20 > nor02_001_sr_qual20.histo &
#${ID}
#30
jellyfish count -m 21 -s 100M -Q ? -t 10 -C <(zcat /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_1.fq.gz) <(zcat /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_2.fq.gz) -o ${ID}_sr_qual30 &
jellyfish histo -t 10 --high=1000000 ${ID}_sr_qual30 > ${ID}_sr_qual30.histo
#20
jellyfish count -m 21 -s 100M -Q 5 -t 10 -C <(zcat /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_1.fq.gz) <(zcat /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_2.fq.gz) -o ${ID}_sr_qual20 &
jellyfish histo -t 10 --high=1000000 ${ID}_sr_qual20 > ${ID}_sr_qual20.histo


scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/jellyfish/${ID}_hifi_jf_reads.histo ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/${ID}/jellyfish/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/nor02_001/jellyfish/nor02_001_hifi_jf_reads.histo
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/nor02_001/short_reads/jellyfish/nor02_001_sr_qual30.histo
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/nor02_001/short_reads/jellyfish/nor02_001_sr_qual20.histo ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/jellyfish/nor02_001_hifi/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/jellyfish_histos.txt ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/jellyfish/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/short_reads/${ID}_sr_qual30.histo  ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/jellyfish/${ID}_hifi/

scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag/ragtag_onlyeight/ragtag_output/SE02_002_corrected_contigs_against_v5.paf ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/SE02_002
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag/ragtag_correction_with_reads_updated_ref/ragtag_without_LR_updated_ref/SE02_scaffolds_witch_correction_LR_updated_ref.paf ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/SE02_002
scp marimond@dell-node-12.mpipz.mpg.de:/biodata/dep_coupland/common/jasmin/male/isabella/stats_09 ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/
# then load histogram in genomescope with k = m, read length = 150, max kmer cov = 1000000
scp  marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/protein_sequences/OrthoFinder/Results_Sep30_1/ ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/Results/


############################################################################################################################################################
#Start Assembly#
ssh marimond@hpcws001.mpipz.mpg.de

#lsfenv
#bjobs
#bqueuesFONT

#${ID}_hifi on hcp cluster

####FLYE####
#go into output directory
flye --version        #2.9-b1768 (19.06.2023)
bsub -q bigmem -R "rusage[mem=200000]" -M 250000 "flye --pacbio-hifi ${col2} -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/flye -g 384564401 -t 45"
bsub -q bigmem -R "rusage[mem=200000]" -M 250000 "flye --pacbio-hifi /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_08_10/m64078e_220808_101749.hifi_reads.fastq.gz -g 420000000 -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/GR01_001/flye -t 50"
bsub -q bigmem -R "rusage[mem=200000]" -M 250000 "flye --pacbio-hifi /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210716_142718.hifi_reads.fastq.gz -g 368043576 -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/nor02_001/pacbio_hifi/flye -t 50"
bsub -q bigmem -R "rusage[mem=200000]" -M 250000 "flye --pacbio-hifi /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz -g 360043576 -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/pacbio_hifi/flye -t 50"
bsub -q bigmem -R "rusage[mem=260000]" -M 300000 "flye --pacbio-hifi /netscratch/dep_coupland/grp_fulgione/danijel/globus/891_2_CELL1/m54274Ue_230623_095607.bc1002--bc1002.hifi_reads.fastq.gz -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/IT19_019/flye -g 384564401 -t 60"

####SMARTDENOVO#####
/netscratch/dep_coupland/grp_fulgione/male/smartdenovo/smartdenovo.pl -t 20 -c 1 -p SE02_002 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz  > SE02_002_hifi_smdn.mak
bsub -q multicore20 -R "rusage[mem=20000]" -M 35000 "make -f SE02_002_hifi_smdn.mak"
/netscratch/dep_coupland/grp_fulgione/male/smartdenovo/smartdenovo.pl -t 20 -c 1 -p NO02_001 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210716_142718.hifi_reads.fastq.gz > NO02_001_hifi_smdn.mak
bsub -q multicore20 -R "rusage[mem=20000]" -M 35000 "make -f NO02_001_hifi_smdn.mak"
####HifiASM####
hifiasm --version     #0.7
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/
bsub -q bigmem -R "rusage[mem=120000]" -M 130000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/LR_Gapcloser/iteration-3/gapclosed.fasta --frags AT03_004_shortread_sort.bam --threads 20 --outdir pilon"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/CH05_002/
bsub -q bigmem -R "rusage[mem=120000]" -M 130000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/CH05_002/LR_Gapcloser/iteration-3/gapclosed.fasta --frags CH05_002_shortread_sort.bam --threads 20 --outdir pilon"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/CZ01_005/
bsub -q bigmem -R "rusage[mem=120000]" -M 130000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/CZ01_005/LR_Gapcloser/iteration-3/gapclosed.fasta --frags CZ01_005_shortread_sort.bam --threads 20 --outdir pilon"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES24_001/
bsub -q bigmem -R "rusage[mem=150000]" -M 160000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES24_001/LR_Gapcloser/iteration-3/gapclosed.fasta --frags ES24_001_shortread_sort.bam --threads 20 --outdir pilon"

bsub -q bigmem -R "rusage[mem=120000]" -M 130000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/CH05_002/LR_Gapcloser/iteration-3/gapclosed.fasta --frags ES03_014_shortread_sort.bam --threads 20 --outdir pilon"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o GR01_001_hifi.asm -t 40 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_08_10/m64078e_220808_101749.hifi_reads.fastq.gz"

#Hifiasm Downstream
#marimond@hpc001
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load python/v3.8.0
conda create --name liftoff


#version 0.16.1
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o GR01_001_hifiasm_0.19.asm -t 40 /netscratch/dep_coupland/grp_fulgione/danijel/globus/891_6_CELL_merged/GR01_001_merged.fq.gz"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o SI01_005_hifiasm_0.19.asm -t 40 /netscratch/dep_coupland/grp_fulgione/danijel/globus/891_1_CELL1/m54274Ue_230705_214109.bc1001--bc1001.hifi_reads.fastq.gz -l0"

seqtk seq -a /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210717_202650.hifi_reads.fastq.gz  > NO01_005_raw.fasta
seqtk seq -a /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210729_191323.hifi_reads.fastq.gz > ES18_040_raw_reads.fasta
seqtk seq -a /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078_210408_170302.Q20.fastq.gz > AT03_004_raw.fasta
seqtk seq -a /netscratch/dep_coupland/grp_fulgione/danijel/globus/891_1_CELL1/m54274Ue_230705_214109.bc1001--bc1001.hifi_reads.fastq.gz > SI01_005_raw.fasta
export PATH="/netscratch/dep_coupland/grp_fulgione/male/assemblies/LR_Gapcloser/src/:$PATH"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "bash LR_Gapcloser.sh -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/SI01_005/RagTag_minimap/ragtag_output/ragtag.scaffold.fasta -l /netscratch/dep_coupland/grp_fulgione/male/assemblies/SI01_005/SI01_005_raw.fasta -m 1000 -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/SI01_005/LR_Gapcloser_minimap"
########ONT Integration in hifiasm########
#assemblies performed with Hifiasm (31) only needed the specification of a parameter for small genomes (-f0)
#and the disabling of purging of duplicated contigs recommended for inbred genomes (-l0)
#Nor2
#wo
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o IT19_017_hifiasm_0.19.asm -t 40  /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078_210408_170302.Q20.fastq.gz"

source /opt/share/software/scs/appStore/modules/init/profile.sh
module load hifiasm/v0.19
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ES24_001_hifiasm_0.19_purged.asm -t 40 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210801_073647.hifi_reads.fastq.gz"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ES20_028_hifiasm_0.19_purged.asm -t 40 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210729_005338.hifi_reads.fastq.gz"
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ES06_006_hifiasm_purged -t 40 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210727_184141.hifi_reads.fastq.gz /biodata/dep_coupland/grp_fulgione/rawdata/alpina_21_12_2020/F20FTSEUHT1304_ARAoxrR/Clean/6/FP100001891BR_L01_SP2012010034-551_1.fq.gz /biodata/dep_coupland/grp_fulgione/rawdata/alpina_21_12_2020/F20FTSEUHT1304_ARAoxrR/Clean/6/FP100001891BR_L01_SP2012010034-551_2.fq.gz"

mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/NO02_001/hifiasm_0.19_ONT/
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/NO02_001/hifiasm_0.19_ONT/
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o NO02_001_hifiasm_0.19_ONT.asm -l0 --ul /biodata/dep_coupland/grp_fulgione/rawdata/alpina_ONT/4959.B/g0103.tar /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210716_142718.hifi_reads.fastq.gz"

mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES03_014/hifiasm_0.19_ONT/
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES03_014/hifiasm_0.19_ONT/
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ES03_014_hifiasm_0.19_ONT.asm -l0 --ul /biodata/dep_coupland/grp_fulgione/rawdata/alpina_ONT/4959.Q/g0097.tar /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210731_012507.hifi_reads.fastq.gz"

mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/hifiasm_0.19_ONT/
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/hifiasm_0.19_ONT/
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o SE02_002_hifiasm_0.19_ONT.asm -l0 --ul /biodata/dep_coupland/grp_fulgione/rawdata/alpina_ONT/4959.F/g0100.tar /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz"

mkdir /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES06_006/hifiasm_0.19_ONT/
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES06_006/hifiasm_0.19/
bsub -q multicore40 -R "rusage[mem=60000]" -M 70000 "hifiasm -o ES06_006_hifiasm_0.19_ONT.asm -l0 --ul /biodata/dep_coupland/grp_fulgione/rawdata/alpina_ONT/4959.S/g0104.tar /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210727_184141.hifi_reads.fastq.gz"

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/hifiasm/${ID}_hifiasm.p_ctg.fa
#produce fasta
awk '/^S/{print ">"$2;print $3}' SE02_002_hifiasm_0.19.asm.bp.p_ctg.gfa >  SE02_002_hifiasm_0.19.asm.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' SE02_002_hifiasm_0.19.asm.bp.hap1.p_ctg.gfa >  SE02_002_hifiasm_0.19.asm.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' SE02_002_hifiasm_0.19.asm.bp.hap2.p_ctg.gfa >  SE02_002_hifiasm_0.19.asm.bp.hap2.p_ctg.fa

#NOR2_001 hifiasm_0.19
awk '/^S/{print ">"$2;print $3}' AT03_004_hifiasm_0.19.asm.bp.hap1.p_ctg.gfa > AT03_004_hifiasm_0.19.asm.bp.hap1.p_ctg.fa &
awk '/^S/{print ">"$2;print $3}' Nor2_01_hifiasm_0.19.asm.bp.hap2.p_ctg.gfa > Nor2_01_hifiasm_0.19.asm.bp.hap2.p_ctg.fa &
awk '/^S/{print ">"$2;print $3}' Nor2_01_hifiasm_0.19.asm.bp.p_ctg.gfa > Nor2_01_hifiasm_0.19.asm.bp.p_ctg.fa &

#NOR2_001 hifiasm_0.19
awk '/^S/{print ">"$2;print $3}' ${col1}_hifiasm_0.19.asm.bp.hap1.p_ctg.gfa  > ${col1}_hifiasm_0.19.asm.bp.hap1.p_ctg.fa  &
awk '/^S/{print ">"$2;print $3}' ${col1}_hifiasm_0.19.asm.bp.hap2.p_ctg.gfa > ${col1}_hifiasm_0.19.asm.bp.hap2.p_ctg.fa &
awk '/^S/{print ">"$2;print $3}' SE06_020_hifiasm_0.19.bp.p_ctg.gfa > SE06_020_hifiasm_0.19.bp.p_ctg.fa
#START QUAST AND BUSCO
https://busco-archive.ezlab.org/frame_plants.html used embryophyta odb10
/opt/share/software/packages/miniconda3-4.12.0/bin/activate bases
#conda install -c bioconda --name bases quast
#Version: 5.0.2

/home/marimond/.conda/envs/bases/bin/quast.py /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.p_ctg.fa -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/hifiasm/quast/danijel -e -f -b
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/hifiasm/${ID}_hifiasm.p_ctg.fa -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/hifiasm_0.19/quast/ -e -f -b &

/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/flye/assembly.fasta -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/flye/quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/flye/assembly.fasta -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/flye/quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/smartdenovo/${ID}_hifi.dmo.cns -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/smartdenovo/quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.19.asm.bp.hap1.p_ctg.fa -o hap1_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.19.asm.bp.hap2.p_ctg.fa -o hap2_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.19.asm.bp.p_ctg.fa -o both_quast -e -f -b &

/home/marimond/.conda/envs/bases/bin/quast.py Nor2_01_hifiasm_0.19.asm.bp.hap1.p_ctg.fa -o hap1_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py Nor2_01_hifiasm_0.19.asm.bp.hap2.p_ctg.fa -o hap2_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py Nor2_01_hifiasm_0.19.asm.bp.p_ctg.fa -o both_quast -e -f -b &

/home/marimond/.conda/envs/bases/bin/quast.py GR01_001_hifi.asm.bp.hap1.p_ctg.fa -o hap1_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py GR01_001_hifi.asm.bp.hap2.p_ctg.fa -o hap2_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py GR01_001_hifi.asm.bp.p_ctg.fa -o both_quast -e -f -b &

/home/marimond/.conda/envs/bases/bin/quast.py SE02_002_hifiasm_0.19.asm.bp.hap1.p_ctg.fa -o hap1_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py SE02_002_hifiasm_0.19.asm.bp.hap2.p_ctg.fa -o hap2_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py  merged_out.fasta -o quast_merged_purged -e -f -b &

busco -i merged_out -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_merged_purged -m genome -f &
busco -i SE02_002_hifiasm_0.19.asm.bp.hap2.p_ctg.fa -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_hap2 -m genome -f &
busco -i ES20_028_final_assembly_minimap.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_ES24_001 -m genome -e

Nor2_01_hifiasm_0.19_ONT.asm.bp.hap1.p_ctg.fa
#BUSCO
#before purging
#Danijel it20 (without purging)
busco -i /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.p_ctg.fa -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats_it20_007 -m genome -f
#Male it20
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/hifiasm/quast/male/
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/hifiasm/it20_007_hifiasm.p_ctg.fa -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome -f

#hifiasm_0.16
#All
/home/marimond/.conda/envs/bases/bin/quast.py smartdn_it20_hapdup_out/hapdup_phased_1.fasta -o hap1_quast -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py merged_GR01_001_hap2_merged.fasta -o quast_hap2_merged -e -f -b &

busco -i merged_GR01_001_hap2_merged.fasta -o busco_hap2_merged -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
#Hap1
busco -i smartdn_it20_hapdup_out/hapdup_phased_2.fasta -o hap2_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
#Hap2
busco -i it20_007_hifiasm_0.19.asm.bp.p_ctg.fa -o both_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &

#Nor2_01_hifiasm_0.19
busco -i GR01_001_hifi.asm.bp.hap1.p_ctg.fa -o hap1_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
#Hap1
busco -i GR01_001_hifi.asm.bp.hap2.p_ctg.fa  -o hap2_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
#Hap2
busco -i GR01_001_hifi.asm.bp.p_ctg.fa  -o both_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &


#flye
cd cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/quast
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/assembly.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome
#SMARTDENOVO
cd cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/quast
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/it20_007_hifi.dmo.cns -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome

#################################################################
#After Purging
#HifiASM
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/purge_haplotigs/quast_purged_hifiasm_results
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/purge_haplotigs/it20_007_purged_flye_pb.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome -f

#FLYE
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/purge_haplotigs/quast_results_purged_flye
#SMARTDENOVO
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/purge_haplotigs
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/purge_haplotigs/it20_007_purged_smartdenovo_pb.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome -f

#after merging
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/quast_hifiasm_flye
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/it20_007_merged_hifiasm_flye_pb.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome -f

/home/marimond/.conda/envs/bases/bin/quast.py mv purged.fa SE02_002_hifiasm_only_purged.fastaa -o quast_both_purged -e -f -b &

busco -i mv purged.fa SE02_002_hifiasm_only_purged.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_both_purged -m genome -f &
###################################################################################################
#after HapDup
#Hifiasm male
#Hap1
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/quast_hifiasm_flye
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/it20_007_merged_hifiasm_flye_pb.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome -f
#Hap2
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/quast_hifiasm_flye
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/it20_007_merged_hifiasm_flye_pb.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats -m genome -f

#FLYE
#Hap1
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/hapdup_out_hifiasm/
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/hapdup_out_hifiasm/hapdup_phased_1.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_hap1 -m genome -f &
#Hap2
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/hapdup_out_hifiasm/
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/hapdup_out_hifiasm/hapdup_phased_2.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats_hap2 -m genome -f &

#Hifiasm Danijel
busco -i it20_007_hifiasm_pb_danijel.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats_hap2 -m genome -f &

##############################################################################################################################################################################################################
#PURGE_HAPLOTIGS
conda install -c bioconda -n bases purge_haplotigs
/home/marimond/.conda/envs/bases/bin/purge_haplotigs test #Purge Haplotigs v1.1.2

#hifiasm
#Male
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 10 -ax map-hifi /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/hifiasm/${ID}_hifiasm.p_ctg.fa ${path_raw_LR_data} --secondary=no \
      | samtools sort -m 1G -o ${ID}_hifiasm_pb_aligned.bam"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 15 -ax asm10 /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly_minimap_purged/${col1}_final_assembly_minimap.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/pangenome/gene_pangenome_library/pangenome_transcripts_library_final.fasta | samtools sort -m 1G -o ${col1}_genes_against_assembly.bam"

#danijel
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 10 -ax map-hifi /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.p_ctg.fa  /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_03_16/demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz --secondary=no \
      | samtools sort -m 1G -o it20_007_danijel_hifiasm_pb_aligned.bam -T it20_007_danijel_hifiasm_pb_tmp.ali"

#flye
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -a ../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta unaligned_contigs.fasta > ref_cont.sam"
      GR01_001_merged.fasta
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi ragtag.scaffolds.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz --secondary=no \
            | samtools sort -m 1G -o SE02_002_Chromosomes_reads.bam"
minimap2 -t 20 -ax map-hifi pilon.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz --secondary=no \
            | samtools sort -m 1G -o SE02_002_pilon_long_reads.bam

            bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi ragtag.scaffold.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz --secondary=no \
                        | samtools sort -m 1G -o SE02_002_Chromosomes_reads_correction_improved.bam"

                        # Using seqtk
seqtk sam2fq ref_cont_sort.sam | seqtk seq -A | \
paste - - | awk '{print $1, $2}' | \
awk 'BEGIN {total=0; match=0} $1==$2 {match++} {total++} END {print "Identity:", (match/total)*100, "%"}'

############################################################################################################################################################
#Purge_dups
#Map reads to assembly and calculate depth +
minimap2 -xasm20 SE02_002_hifiasm_0.19.asm.bp.p_ctg.fa /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz   | gzip -c - > SE02_002_hifiasm_pb.paf.gz
pbcstat GR01_001_both_hifiasm_pb.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

#split assembly and do self-alignment
split_fa GR01_001_hifi.asm.bp.p_ctg.fa > GR01_001_both.split
minimap2 -xasm5 -DP  GR01_001_both.split GR01_001_both.split | gzip -c - >  GR01_001_both.split.self.paf.gz


#purge haplotigs and overlaps
purge_dups -2 -T cutoffs -c GR01_001_both.split.self.paf.gz > dups.bed 2> purge_dups.log

#get purged sequences and haplotigs
get_seqs dups.bed GR01_001_hifi.asm.bp.p_ctg.fa


bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bash purge_dups.sh"
#creating a coverage histogram
#hifiasm male
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b it20_007_hifiasm_0.16_hap1_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.16/it20_007_hifiasm_0.16.asm.bp.hap1.p_ctg.fa -t 10 &
#hifiasm Danijel
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b it20_007_danijel_hifiasm_pb_aligned.bam -g /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.p_ctg.fa -t 10 &
#smartdenovo
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist -b /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/purge_haplotigs/it20_007_smartdenovo_pb_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/it20_007_hifi_new.dmo.cns -t 10
#flye
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist -b /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/purge_haplotigs/it20_007_flye_pb_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/assembly.fasta -t 10

/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b it20_007_hifiasm_0.19_hap2_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.16/it20_007_hifiasm_0.16.asm.bp.hap2.p_ctg.fa -t 10 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b it20_007_hifiasm_0.19_both_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.16/it20_007_hifiasm_0.16.asm.bp.both.p_ctg.fa -t 10 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b it20_007_hifiasm_0.19_hap1_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.16/it20_007_hifiasm_0.16.asm.bp.hap1.p_ctg.fa -t 10 &
#flye after hapDup
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b it20_007_flye_hap1_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/hapdup_out_flye/hapdup_phased_1.fasta -t 10 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b it20_007_flye_hap2_aligned.bam -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/hapdup_out_flye/hapdup_phased_2.fasta -t 10 &

/home/marimond/.conda/envs/bases/bin/purge_haplotigs  hist  -b SE02_002_flye_pb_aligned.bam -g assembly.fasta -t 10 &
n9T9dJ7ppn
#Manual: Set coverage limits
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.16/it20_007_hifiasm_0.16_both_aligned.bam.histogram.png ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/it20_007_hifi/hifiasm_0.16/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/purge_haplotigs/it20_007_flye_pb_aligned.bam.histogram.png ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/it20_007_hifi/flye/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/purge_haplotigs/it20_007_hifiasm_pb_aligned.bam.histogram.png ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/it20_007_hifi/hifiasm/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/purge_haplotigs/it20_007_smartdenovo_pb_aligned.bam ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/it20_007_hifi/smartdenovo/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/purge_haplotigs/danijel/it20_007_danijel_hifiasm_pb_aligned.bam.histogram.png ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/it20_007_hifi/flye/

scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.19/it20_007_hifiasm_0.19_both_aligned.bam.histogram.png ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/it20_007_hifi/flye/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/hapdup_out_flye/it20_007_flye_hap2_aligned.bam.histogram.png ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/it20_007_hifi/flye/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/nor02_001/pacbio_hifi/flye/nor02_001_flye_pb_aligned.bam.histogram.png ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/nor02_001/
scp marimond@hpc001:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag/ragtag_output_old/ragtag_output/SE02_scaffolds_v5.1.bam ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/SE02_002

#Run Purge_haplotigs coverage statistics
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_smartdenovo_pb_aligned.bam.gencov  -l 6  -m 27 -h 70  -o it20_007_hifiasm_pb_coverage_stats.csv -j 80 -s 80
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_flye_pb_aligned.bam.gencov  -l 5  -m 32 -h 60  -o it20_007_flye_pb_coverage_stats.csv -j 80 -s 80
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_hifiasm_pb_aligned.bam.gencov  -l 5  -m 26 -h 75  -o it20_007_hifiasm_pb_coverage_stats.csv -j 80 -s 80
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_danijel_hifiasm_pb_aligned.bam.gencov  -l 5  -m 22 -h 75  -o it20_007_hifiasm_pb_coverage_stats_danijel.csv -j 80 -s 80
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_hifiasm_0.16_hap2_aligned.bam.gencov  -l 5  -m 23 -h 175  -o it20_007_hifiasm_0.16_hap2.csv -j 80 -s 80 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_hifiasm_0.16_hap1_aligned.bam.gencov  -l 5  -m 25 -h 175  -o it20_007_hifiasm_0.16_hap1.csv -j 80 -s 80
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_hifiasm_0.16_both_aligned.bam.gencov  -l 5  -m 25 -h 175  -o it20_007_hifiasm_0.16_both.csv -j 80 -s 80

/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_hifiasm_0.19_both_aligned.bam.gencov  -l 5  -m 25 -h 175  -o it20_007_hifiasm_0.19_both.csv -j 80 -s 80 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_hifiasm_0.19_hap1_aligned.bam.gencov -l 5  -m 24 -h 175  -o it20_007_hifiasm_0.19_hap1.csv -j 80 -s 80 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_hifiasm_0.19_hap2_aligned.bam.gencov -l 4  -m 24 -h 175  -o it20_007_hifiasm_0.19_hap2.csv -j 80 -s 80 &
#FLYE
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_flye_hap2_aligned.bam.gencov-l 2  -m 33 -h 180  -o it20_007_flye_hap2.csv -j 80 -s 80 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs cov  -i it20_007_flye_hap1_aligned.bam.gencov -l 2  -m 33 -h 180  -o it20_007_flye_hap1.csv -j 80 -s 80 &



#Finally run the purging procedure
purge_haplotigs  purge  -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/hifiasm/${ID}_hifiasm.p_ctg.fa  -c ${ID}_hifiasm_pb_coverage_stats.csv -o ${ID}_hifiasm_pb
purge_haplotigs  purge  -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/hifiasm/${ID}_hifiasm.p_ctg.fa  -c ${ID}_hifiasm_pb_coverage_stats.csv -o ${ID}_hifiasm_pb
purge_haplotigs  purge  -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/${ID}/pacbio_hifi/hifiasm/${ID}_hifiasm.p_ctg.fa  -c ${ID}_hifiasm_pb_coverage_stats.csv -o ${ID}_hifiasm_pb

purge_haplotigs purge -g it20_007_hifiasm_0.19.asm.bp.p_ctg.fa -c it20_007_hifiasm_0.19_both.csv -o it20_007_hifiasm_0.19_both_purged &
purge_haplotigs purge -g it20_007_hifiasm_0.19.asm.bp.hap1.p_ctg.fa -c it20_007_hifiasm_0.19_hap1.csv -o it20_007_hifiasm_0.19_hap1_purged &
purge_haplotigs purge -g it20_007_hifiasm_0.19.asm.bp.hap2.p_ctg.fa -c it20_007_hifiasm_0.19_hap2.csv -o it20_007_hifiasm_0.19_hap2_purged &

purge_haplotigs purge -g hapdup_phased_2.fasta -c it20_007_flye_hap2.csv -o it20_007_flye_hap2_purged &
purge_haplotigs purge -g hapdup_phased_1.fasta -c it20_007_flye_hap1.csv -o it20_007_flye_hap1_purged






/home/marimond/.conda/envs/bases/bin/purge_haplotigs purge  -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/assembly.fasta  -c it20_007_flye_pb_coverage_stats.csv -o it20_007_purged_flye_pb &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs purge  -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/it20_007_hifi_new.dmo.cns  -c it20_007_hifiasm_pb_coverage_stats.csv -o it20_007_purged_smartdenovo_pb
/home/marimond/.conda/envs/bases/bin/purge_haplotigs purge  -g /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/it20_007_hifiasm.p_ctg.fa  -c it20_007_hifiasm_pb_coverage_stats.csv -o it20_007_hifiasm_pb
/home/marimond/.conda/envs/bases/bin/purge_haplotigs purge  -g /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.p_ctg.fa  -c it20_007_hifiasm_pb_coverage_stats_danijel.csv -o it20_007_hifiasm_pb_danijel

/home/marimond/.conda/envs/bases/bin/purge_haplotigs purge  -g it20_007_hifiasm_0.16.asm.bp.hap2.p_ctg.fa  -c it20_007_hifiasm_0.16_hap2.csv -o it20_007_hifiasm_0.16_purged_hap2 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs purge  -g it20_007_hifiasm_0.16.asm.bp.hap1.p_ctg.fa -c it20_007_hifiasm_0.16_hap1.csv -o it20_007_hifiasm_0.16_purged_hap1 &
/home/marimond/.conda/envs/bases/bin/purge_haplotigs purge  -g it20_007_hifiasm_0.16.asm.bp.p_ctg.fa -c it20_007_hifiasm_0.16_both.csv -o it20_007_hifiasm_0.16_purged_both &


busco -i assembly.fasta -o both_purged_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb1 -m genome -f &
#Hap1
busco -i it20_007_hifiasm_0.16_purged_hap1.fasta -o hap1_purged_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
#Hap2
busco -i it20_007_hifiasm_0.16_purged_hap2.fasta -o hap2_purged_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &

######WORK WITH PHASED HAPS?############
redundans.py -f ../genome.fa -t 32 --nocleaning --noscaffolding \
 --norearrangements --nogapclosing

#Quality assesment of purged sets with
#QUAST
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/purge_haplotigs/it20_007_hifiasm_pb.fasta -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/purge_haplotigs/it20_007_flye_pb.fasta -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/purge_haplotigs/it20_007_purged_smartdenovo_pb.fasta -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/purge_haplotigs/danijel/it20_007_hifiasm_pb_danijel.fasta -e -f -b &

/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.19_hap1_purged.fasta -o quast_hap1_purged -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.19_hap2_purged.fasta -o quast_hap2_purged -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.19_both_purged.fasta -o quast_both_purged -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py hapdup_phased_2.fasta -o quast_hap2_flye -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py hapdup_phased_1.fasta -o quast_hap1_flye -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py Nor02_ONT_polished_sr.fasta -o quast_ONT_nor02_polished -e -f -b &
busco -i Nor02_ONT_polished_sr.fasta -o busco_ONT_nor02_polished -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
busco -i it20_007_hifiasm_0.19_hap2_purged.fasta -o hap2_purged_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
busco -i it20_007_hifiasm_0.19_both_purged.fasta -o both_purged_busco -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
busco -i hapdup_phased_1.fasta -o busco_hap1_flye -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &
busco -i hapdup_phased_2.fasta -o busco_hap2_flye -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -m genome -f &




#################################################################################################################################
#PHASE ASSEMBLIES WITH HAPDUP#
minimap2 -t 20 -ax map-hifi assembly.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_08_10/m64078e_220808_101749.hifi_reads.fastq.gz --secondary=no \
      | samtools sort -m 1G -o GR01_001_flye_pb_aligned.bam


singularity pull docker://mkolmogo/hapdup:0.12
HD_DIR=`pwd`

#flye
singularity exec --bind $HD_DIR hapdup_0.12.sif \
  hapdup --assembly assembly.fasta \
   --bam GR01_001_flye_pb_aligned.bam \
   --bam-index GR01_001_flye_pb_aligned.bam.bai \
   --out-dir GR01_001_hapdup_out -t 10 --rtype hifi
#hifiasm
singularity exec --bind $HD_DIR hapdup_0.12.sif \
  hapdup --assembly $HD_DIR/it20_007_hifi.dmo.cns.fasta \
   --bam $HD_DIR/it20_007_smartdenovo_pbreads_aligned.bam \
   --bam-index $HD_DIR/it20_007_smartdenovo_pbreads_aligned.bam.bai \
   --out-dir $HD_DIR/smartdn_it20_hapdup_out -t 10 --rtype hifi

#QUAST
#Quality assesment of phased sets
#flye
/home/marimond/.conda/envs/bases/bin/quast.py ../hapdup_phased_1.fasta -o quast_hapdup_phased_1 -e -f -b
/home/marimond/.conda/envs/bases/bin/quast.py ../hapdup_phased_2.fasta -o quast_hapdup_phased_2 -e -f -b &
#danijel
/home/marimond/.conda/envs/bases/bin/quast.py /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.hap1.p_ctg.fa -o danijel_hap1 -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.hap2.p_ctg.fa -o danijel_hap2 -e -f -b &
#hifiasm Male
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/hapdup_out_hifiasm/hapdup_phased_1.fasta -o hap1 -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/hapdup_out_hifiasm/hapdup_phased_2.fasta -o hap2 -e -f -b &

#BUSCO
#Danijels hpa1/2 hifiasm
cd
busco -i /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.hap1.p_ctg.fa -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats_hap1 -m genome -f
busco -i /biodata/dep_coupland/grp_fulgione/danijel/assemblies/self_incompatible_italians/ital_11_5440C_IT2-07/ital_3.bp.hap2.p_ctg.fa -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats_hap2 -m genome -f
#flye
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/hapdup_out_flye/hapdup_phased_1.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats_hap1 -m genome -f &
busco -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/hapdup_out_flye/hapdup_phased_2.fasta -l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_stats_hap2 -m genome -f &

##################################################################################################################################
#Merge purged assembly sets with quickmerge
#nucmer 4.0.0beta2 # Nucmer aligns the two assemblies so that the merger can find the correct splice sites:
# assembly with the longest contiguity (N50) selected as query, the assembly with second longest N50 was used as reference to join query assembly
#The resulting assembly was further improved by using the third assembly as reference


#####run  EXPORT#########
export PATH=/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge/MUMmer3.23:$PATH
merge_wrapper.py it20_007_purged_flye_pb.fasta it20_007_hifiasm_0.16.asm.bp.p_ctg.fa

export PATH=/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge/MUMmer3.23:$PATH
bsub -q multicore40 -R "rusage[mem=30000]" -M 40000 "merge_wrapper.py /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES20_028/flye/purged.fa /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES20_028/hifiasm_0.19_purged/ES20_028_hifiasm_0.19.asm.bp.p_ctg.fa -l 2000000 -t 20"
bsub -q multicore40 -R "rusage[mem=30000]" -M 40000 "augustus --species=arabidopsis unaligned_contigs_renamed.fasta > unaln_contigs_annotated"
bsub -q multicore40 -R "rusage[mem=60000]" -M 65000 "cd-hit-para.pl -i unaligned_contigs_renamed.fasta -o unaligned_contigs_non-redundant --P cd-hit-est --L 40"

/netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/flye/purged.fa /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/hifiasm_0.19_purged/AT03_004_hifiasm_0.19.asm.bp.p_ctg.fa -l 200000
paste("/netscratch/dep_coupland/grp_fulgione/male/assemblies/",name,"final_assembly")
/netscratch/dep_coupland/grp_fulgione/male/assemblies/ DE01_001 /final_assembly//chr_chr1_coverage_mean.csv

/home/marimond/.conda/envs/bases/bin/quast.py merged_it20_flye_hifiasm.19.fasta -o quast_merged_new -e -f -b &
busco -i merged_it20_flye_hifiasm.19.fasta-l /netscratch/dep_coupland/grp_fulgione/danijel/busco/brassicales_odb10 -o busco_merged_new -m genome -f &
#quast for merged sets
/home/marimond/.conda/envs/bases/bin/quast.py merged_SE02_002_merge_before.fasta -o quast_merged_wo_purged_se02 -e -f -b &
#Quast for purged sets
/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.16_purged_hap1.fasta -o purged_quast_hap1 -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py it20_007_hifiasm_0.16_purged_hap2.fasta -o purged_quast_hap2 -e -f -b &
/home/marimond/.conda/envs/bases/bin/quast.py merged.fasta -o merged_quast_both -e -f -b &
#merge hap fyle_hifiasm
/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge/merge_wrapper.py

#Now add smartdenovo to the merged hifiasm_flye set (Smartdenovo as reference)
nucmer -p it20_007_hifiasm_flye_smartdenovo /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/purge_haplotigs/it20_007_purged_smartdenovo_pb.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/it20_007_merged_hifiasm_flye_pb.fasta &
delta-filter -r -q -l 10000 it20_007_hifiasm_flye_smartdenovo.delta > it20_007_hifiasm_flye_smartdenovo_rq.delta &
/quickmerge -d it20_007_hifiasm_flye_smartdenovo_rq.delta -r /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/purge_haplotigs/it20_007_purged_smartdenovo_pb.fasta  -q /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merged_assembly/it20_007_merged_hifiasm_flye_pb.fasta -hco 5.0 -c 1.5 -l 1600000 -ml 9000 -p it20_007_hifiasm_flye_smartdenovo

#failes with SMARTDENOVO



######################################################################################################################################
#GALA Gap free
git clone https://github.com/ganlab/gala.git
cd GALA
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "gala/gala -a hifiasm flye smartdenovo -f out_gala draft_names_paths.txt fq /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_03_16/demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz pacbio-raw"
cp gala/src/read_extract.py gala/src/save_original_read_extract.py
awk '{gsub("\t","    "); print}' gala/src/save_original_read_extract.py > gala/src/read_extract.py
#1) prepare draft txt
####GALA Running as gala/gala input -pabio-raw -> on cluster
touch draft_names_paths.txt
draft_01=/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm/it20_007_hifiasm.p_ctg.fa
draft_02=/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/flye/assembly.fasta
draft_03=/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/smartdenovo/it20_007_hifi_new.dmo.cns

#MISASSEMBLY MODULE

#1) Use the comp module to generate a draft_comparison file
gala/comp draft_names_paths.txt

#2) Run draft_comparison file to produce drafts comparison paf files
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "sh draft_compare.sh"

#3) Use the mdm module to identify mis-assembled contigs.
gala/mdm comparison/ 3

#4) Use the newgenome module to Produce misassembly-free drafts.
gala/newgenome draft_names_paths.txt gathering/

#CONTIG CLUSTERING MODULE
#5) Use the comp module to generate a draft_comparison file for misassembly-free drafts.
gala/comp new_draft_names_paths.txt

#6) Run draft_comparison file to produce new drafts comparison paf files.
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bash draft_comp.sh"

#7) Run the ccm module to produce contigs scaffolding groups.
gala/ccm comparison_new/ 3


#SCAFFOLDING GROUP ASSEMBLY
#8) Map all drafts against raw long reads and self-corrected reads if available.
bwa index new_draft_01.fa
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bwa mem -x pacbio new_draft_01.fa /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_03_16/demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz | samtools sort | samtools view -Sb > new_draft_01_mapping.bam"
bwa index new_draft_02.fa
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bwa mem -x pacbio new_draft_01.fa /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_03_16/demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz | samtools sort | samtools view -Sb > new_draft_02_mapping.bam"
bsub -q multicore20 -R "rusage[mem=25000]" -M 30000 "index new_draft_03.fa"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "bwa mem -x pacbio new_draft_01.fa /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/HiFi_2022_03_16/demultiplex.bc1011_BAK8A_OA--bc1011_BAK8A_OA.hifi_reads.fastq.gz | samtools sort | samtools view -Sb > new_draft_03_mapping.bam"


#9) Use the following commands to separate the read names mapped to each contig

samtools view -H new_draft_01_mapping.bam |grep "SQ"|cut -f 2|cut -d : -f 2 > 01_new_draft_contig_names

seprator 01_new_draft_contig_names new_draft_01_mapping.bam

sh bam_seprator.sh

for i in bams/*; do samtools view $i | cut -f 1 > $i.read_names;done;




All the results are only for one genome (it20_007)
-so only of the comparison of the generated hifiasm set by Danijel and me it seems that after purging the differences in genome length between the two sets gets less. But the N50 of the newer hifiasm set is still a bit higher.
-After phasing my Hifiasm with HapDup however the total length and N50 and BUSCO of Danijel default phased Haplotypes show better results. So if we decide to go with the phased haplotype for the outcrossers I would suggest that we go with the hifiasm haps Danijel created?:)
-Flye and hifiasm show both good results in BUSCO but flye produces much more contigs also after purging then hifiasm, when I am back I wanna try the scaffolding on both and see how it changes.
-Smartdenovo shows form BUSCO as well as from quast not that good results compared tp the other, so maybe we shouldnt include it in the merging set?

-GALA is still running so I wait for the results



############################################################################################
#Scaffolding# with RagTag
export PATH="/home/marimond/.conda/envs/bases/bin/:$PATH"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta ../merged_flye_hifiasm_purged/merged_out.fasta -o ragtag_wo_corr -C --aligner /home/marimond/.conda/envs/bases/bin/nucmer -t 20"
# correct a query assembly
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta ../merged_flye_hifiasm_purged/merged_out.fasta -o ES24_001_ragtag_purged -C -t 20"


export PATH="/home/marimond/.conda/envs/bases/bin/:$PATH"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/SI01_005/merged_flye_hifiasm/merged_out.fasta -C -t 20"


bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py correct ../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta ../merged_flye_hifiasm/merged_out.fasta -R /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210722_002416.hifi_reads.fastq.gz -T corr --aligner /home/marimond/.conda/envs/bases/bin/nucmer -o ragtag_nucmer_correction -t 20"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/merged_flye_hifiasm/merged_out.fasta -o ragtag_wo_corr -C --aligner /home/marimond/.conda/envs/bases/bin/nucmer -t 20"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py scaffold ../../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta ragtag.correct.fasta -C -o ragtag_scaffold_nucmer --aligner /home/marimond/.conda/envs/bases/bin/nucmer -t 20"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "ragtag.py patch SE02_002_chr1-8.fasta m64078e_210719_121246.hifi_reads.fastq --fill-only -t 20 -o patch_only_fill"
# make joins and fill gaps in target.fa using sequences from query.fa


samtools faidx -f
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "/home/marimond/.conda/envs/bases/bin/ragtag.py patch ragtag.scaffolds.fasta m64078e_210719_121246.hifi_reads.fastq -o patch"

git clone https://github.com/tpoorten/dotPlotly.git
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "/opt/share/software/packages/MaSuRCA-3.3.4/bin/chromosome_scaffolder.sh -r /biodata/dep_coupland/common/Arabis_alpina_resource/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta  -q ../merge_flye_hifiasm/pilon/nor02_merged_pilon.fasta -s /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64253e_210716_142718.hifi_reads.fastq.gz -i 95 -t 20"
./dotPlotly/pafCoordsDotPlotly.R \
   -i ../pajares_hifiasm19-cont_it20_007.paf \
   -o ../pajares_hifiasm19-cont_it20_007.plot \
   -m 2000 \
   -q 500000 \
   -k 10 \
   -s -t -l -p 12

####################################################
#Cleaning with PILON
#Self-alignment
# option "--cs" is recommended as paftools.js may need it
minimap2 -cx asm5  > nor02_001_ONT_self_alignment.paf

#Short read alignment with BWA mem
#!/bin/bash
assembly="/netscratch/dep_coupland/grp_fulgione/male/assemblies/CH05_002/LR_Gapcloser/iteration-3/gapclosed.fasta"
sr_fw="/biodata/dep_coupland/grp_fulgione/rawdata/alpina_08_06_2021/F20FTSEUHT1304_ARApoigR/Clean/433/FP100002620TL_L01_SP2103181767-81_1.fq.gz"
sr_rv="/biodata/dep_coupland/grp_fulgione/rawdata/alpina_08_06_2021/F20FTSEUHT1304_ARApoigR/Clean/433/FP100002620TL_L01_SP2103181767-81_2.fq.gz"
col_1="CH05_002"
bwa index ${assembly}

echo "start shortread mapping"
bwa mem ${assembly} ${sr_fw} ${sr_rv} > ${assembly}_shortread.sam > ${col1}_shortread.sam

echo "create bam and index"
samtools view -S -b ${col1}_shortread.sam > ${col1}_shortread.bam
samtools sort ${col1}_shortread.bam -o ${col1}_shortread_sort.bam
samtools index -b  ${col1}_shortread_sort.bam


echo "run pilon"
pilon --genome ${assembly} --frags  ${col1}_shortread_sort.bam --outdir pilon
mv  ${col1}_shortread.bam pilon/
mv ${col1}_shortread_sort.bam pilon/
mv  ${col1}_shortread_sort.bam pilon/
mv  ${col1}_shortread.sam pilon/




#SE02_002
#!/bin/bash
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
output_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/pilon.sh"

while read -r col1 col2 _; do
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}
  echo "
    #!/bin/bash
    bwa index /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser/iteration-3/gapclosed.fasta

    echo "start shortread mapping"
    bwa mem /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser/iteration-3/gapclosed.fasta ${col3} ${col4} > ${col1}_shortread.sam

    echo "create bam and index"
    samtools view -S -b ${col1}_shortread.sam > ${col1}_shortread.bam
    samtools sort ${col1}_shortread.bam -o ${col1}_shortread_sort.bam
    samtools index -b  ${col1}_shortread_sort.bam

    echo "run pilon"
    pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser/iteration-3/gapclosed.fasta --frags ${col1}_shortread_sort.bam --outdir pilon " >> ${output_file}
    mv ${output_file} /netscratch/dep_coupland/grp_fulgione/male/assemblies//${col1/${col1}_pilon.sh
    echo "${col1}_pilon.sh submitted"
    bsub -q multicore40 -R "rusage[mem=65000]" -M 70000 "bash ${col1}_pilon.sh"
done < "$input_file"
echo "all jobs submitted"




#!/bin/bash
col1="ES20_028"
col3="/biodata/dep_coupland/grp_fulgione/rawdata/alpina_08_06_2021/F20FTSEUHT1304_ARApoigR/Clean/669/FP100002644BL_L01_SP2103181964-86_1.fq.gz"
col4="/biodata/dep_coupland/grp_fulgione/rawdata/alpina_08_06_2021/F20FTSEUHT1304_ARApoigR/Clean/669/FP100002644BL_L01_SP2103181964-86_2.fq.gz"

bwa index /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser_minimap_purged/iteration-3/gapclosed.fasta

echo "start shortread mapping"
bwa mem /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser_minimap_purged/iteration-3/gapclosed.fasta ${col3} ${col4} | samtools view -S -b - > ${col1}_shortread.bam
samtools sort ${col1}_shortread.bam -o ${col1}_shortread_sort.bam
samtools index -b ${col1}_shortread_sort.bam

echo "run pilon"

pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/LR_Gapcloser_minimap_purged/iteration-3/gapclosed.fasta --frags ${col1}_shortread_sort.bam --outdir pilon

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT/
bsub -q bigmem -R "rusage[mem=200000]" -M 220000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES18_024/LR_Gapcloser_minimap/iteration-3/gapclosed.fasta --frags ES27_024_shortread_sort.bam --outdir pilon"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/
bsub -q bigmem -R "rusage[mem=200000]" -M 220000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES04_014/LR_Gapcloser_minimap/iteration-3/gapclosed.fasta --frags ES04_014_shortread_sort.bam --outdir pilon"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/
bsub -q bigmem -R "rusage[mem=200000]" -M 220000 "pilon --genome /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/LR_Gapcloser_minimap/iteration-3/gapclosed.fasta --frags AT03_004_shortread_sort.bam --outdir pilon"

cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES20_028/
bsub -q bigmem -R "rusage[mem=200000]" -M 220000 "bash ES20_028_pilon_minimap.sh"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES04_014/
bsub -q bigmem -R "rusage[mem=200000]" -M 220000 "bash ES04_014_pilon_minimap.sh"
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/
bsub -q bigmem -R "rusage[mem=200000]" -M 220000 "bash AT03_004_pilon_minimap.sh"

###############################################################################################
#Self-Alignment with MUMmer
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load 4.0.0rc1

nucmer -p Nor02_merged_ONT_flye -t 16 $ref $query
nucmer -p Nor02_merged_ONT_flye -t 16 /biodata/dep_coupland/common/Arabis_alpina_resource/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta nor02_merged_flye_hifi.fasta
mummerplot -p Nor02_merged_ONT_flye -t png Nor02_merged_ONT_flye.delta



delta-filter -i $IDENTITY -l $SIZE ${PREFIX}.delta > ${PREFIX}_filtered.delta

mummerplot -p ${PREFIX}_filtered -t png ${PREFIX}_filtered.delta



show-coords -r -T -H -d ${PREFIX}.delta > ${PREFIX}.coords





###############G###############


##Post-Analysis of scaffolds
#create paf for genome synteny analysis between reference and sample
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5  ../../../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta ragtag.scaffold.fasta > DE01_scaffolds_with_corr_updated_ref.paf"

bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5  ../../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta ragtag.scaffold.fasta > NO02_scaffolds_wo_correction_updated_ref.paf"
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5  ../../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta ragtag.scaffold.fasta > SE02_scaffolds_wo_correction_updated_ref.paf"
#create bam for IGV analysis between genome and raw reads
bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -a asm5 ../../Reference_v5.1_analysis/Reference_v5.1_updated.fasta unaligned_contigs.fasta --secondary=no \ | samtools sort -m 1G -o unaln_cont_against_ref.bam"
bsub -q multicore40 -R "rusage[mem=65000]" -M 70000 "gapless.sh -j 30 -i ragtag.scaffold.fasta -o gapless -t pb_hifi /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210722_002416.hifi_reads.fastq.gz"
gapless.sh -j 30 -i ragtag.scaffold.fasta -o gapless -t pb_hifi /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210722_002416.hifi_reads.fastq.gz

#For it20_007_purged_smartdenovo_pb
#for IGV
#Read as genome file:
biodata/dep_coupland/common/Arabis_alpina_resource/Arabis_alpina.MPIPZ.version_5.1.chr.all.fasta
#Read in as normal file
biodata/dep_coupland/common/Arabis_alpina_resource/annotation/Arabis_alpina_mpipz_v5.1_annotation_eggnog_genes.gff
netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag/ragtag_output_old/ragtag_output/SE02_002_Chromosomes_reads.bam.bai
netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag/ragtag_output_old/ragtag_output/SE02_002_Chromosomes_reads.bam




#############ILRA POST-Processing of contigs################
/opt/share/software/scs/appStore/stretchApps/buildUtils/mambaforge/v23.1.0-4/condabin/ILRA_env/bin/
#BiocManager::install("BiocGenerics")
#BiocManager::install("BiocVersion")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomeInfoDbData", force=TRUE)
#BiocManager::install("ggtree", force=TRUE)
#BiocManager::install("IRanges", force=TRUE)
#BiocManager::install("S4Vectors", force=TRUE)
#BiocManager::install("treeio", force=TRUE)
#BiocManager::install("XVector", force=TRUE)
#BiocManager::install("zlibbioc", force=TRUE)


#Scaffold Gap closing
#TGSGAPCLOSER
bsub -q bigmem -R "rusage[mem=400000]" -M 450000  "tgsgapcloser  \
        --scaff  SE02_002_chr1-8.fasta \
        --reads  /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078e_210719_121246.hifi_reads.fastq.gz \
        --output tgsgapcloser \
        --racon  /home/marimond/.conda/envs/bases/bin/racon \
        --tgstype pb \
        --threads 40 \
        >pipe.log 2>pipe.err"



  #### Quast evaluation of assemblies ####
/netscratch/dep_coupland/grp_fulgione/male/assemblies/quast-5.2.0/quast.py ragtag.scaffold.fasta --output-dir quast -r /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta

########## GAPCLOSING ###################

/netscratch/dep_coupland/grp_fulgione/male/assemblies/LR_Gapcloser/src$ bash LR_Gapcloser.sh -i ../../SE02_002/RagTag/ragtag_without_correction_updated_ref/ragtag.scaffold.fasta -l ../../SE02_002/SE02_raw.fasta -o ../../LR_Gapcloser/
bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "bash LR_Gapcloser.sh -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/NO02_001/RagTag_nucmer/ragtag_wo_corr/ragtag.scaffold.fasta -l NO02_001_raw.fasta -o LR_Gapcloser/"

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
bjobs


#CREATE PAF AND BAM for Evaluation

#!/bin/bash
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

while read -r col1 col2 _; do
  echo "Submit ${col1}"
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly/
  bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -ax map-hifi AT03_004_final_assembly_minimap.fasta /biodata/dep_coupland/grp_fulgione/rawdata/alpina_Hifi/m64078_210408_170302.Q20.fastq.gz --secondary=no \ | samtools sort -m 1G -o AT03_004_reads_against_assembly.bam"
  bsub -q multicore20 -R "rusage[mem=35000]" -M 45000 "minimap2 -t 20 -cx asm5 /netscratch/dep_coupland/grp_fulgione/male/assemblies/Reference_v5.1_analysis/Reference_v5.1_updated.fasta ${col1}_assembly_only_chr.fasta > ${col1}_scaffolds_against_ref.paf"

done < "$input_file"

echo "all jobs sended"

export PATH="/netscratch/dep_coupland/grp_fulgione/male/assemblies/LR_Gapcloser/src/:$PATH"
bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "bash LR_Gapcloser.sh -i /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/RagTag_nucmer/ragtag_output/ragtag.scaffold.fasta -l /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/AT03_004_raw.fasta -m 1000 -o /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/LR_Gapcloser"

###########Stats##################
input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
samtools_index="/netscratch/dep_coupland/grp_fulgione/male/assemblies/Assembly_stats.txt"
while read -r col1 col2 _; do
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/
  echo "${col1}" >> ${samtools_index}
  echo "Mean Confidence Values"
  cut -f2,3,4 ragtag.scaffold.confidence.txt | awk '{ total_grouping += $1; total_location += $2; total_orientation += $3; count++ } END { print "Mean Grouping Confidence:", total_grouping/count; print "Mean Location Confidence:", total_location/count; print "Mean Orientation Confidence:", total_orientation/count }'  >> ${samtools_index}
   awk /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/ragtag.scaffold.fasta.fai"  >> ${samtools_index}

cut -f2,3,4 ragtag.scaffold.fasta.fai


/biodata/dep_coupland/grp_fulgione/rawdata/alpina_24_03_2021/F20FTSEUHT1304_ARAsovR/Clean/233
#Map short reads#####
bowtie2-build ragtag.scaffold.fasta SE02_bowtie
bowtie2 -x SE02_bowtie -1 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_24_03_2021/F20FTSEUHT1304_ARAsovR/Clean/293/FP100002258TL_L01_SP2101080252-505_1.fq.gz -2 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_24_03_2021/F20FTSEUHT1304_ARAsovR/Clean/293/FP100002258TL_L01_SP2101080252-505_2.fq.gz | samtools view -@ 10 -Sb -o SE_bowtie2.bam
samtools depth -a ${col1}_reads_against_assembly.bam > .depth
bedtools makewindows -g <(cut -f 1,2 gapclosed.fasta.fai) -w 5000 > SE02_002_assembly.bed

/netscratch/dep_coupland/grp_fulgione/male/assemblies/Skripte/RScripts/coverage_analysis_longread.R

#SHORT READ MAPPING
#build index
bowtie2-build /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.19/it20_007_hifiasm_0.19.asm.bp.hap1.p_ctg.fa IT20_007_hap1_bowtie
bowtie2-build /netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/hifiasm_0.19/it20_007_hifiasm_0.19.asm.bp.hap2.p_ctg.fa IT20_007_hap2_bowtie

#map short reads against ccontigs
bowtie2 --threads 20 -x it20_hap1_bowtie -1 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_1.fq.gz -2 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_2.fq.gz | samtools view -@ 10 -Sb -o IT20_007_hap1_short_reads.bam
bowtie2 --threads 20 -x IT20_007_hap2_bowtie -1 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_1.fq.gz -2 /biodata/dep_coupland/grp_fulgione/rawdata/alpina_15_02_2023/clean/865/V350133608_L02_B5GPLAytilRABBA-574_2.fq.gz | samtools view -@ 10 -Sb -o IT20_007_hap2_short_reads.bam

bsub -q multicore40 -R "rusage[mem=85000]" -M 95000 "bash short_read_bowtie_mapping.sh"

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"

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


awk '$1 == "chr1_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr1.coverage
awk '$1 == "chr2_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr2.coverage
awk '$1 == "chr3_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr3.coverage
awk '$1 == "chr4_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr4.coverage
awk '$1 == "chr5_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr5.coverage
awk '$1 == "chr6_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr6.coverage
awk '$1 == "chr7_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr7.coverage
awk '$1 == "chr8_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_chr8.coverage
awk '$1 == "Chloroplast_Aa_NC_023367.1_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_cloroplast.coverage
awk '$1 == "Mitochondrium_Aa_NC_037070.1_RagTag_pilon_pilon" {print $0}' SE02_002_pilon_long_reads.depth > SE02_002_SR_mitochondrium.coverage

scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/pilon/pilon_2/Chromosome1-4_depth_lr_pilon.pdf  ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/SE02_002/final_assembly/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/pilon/pilon_2/Chromosome5-8_depth_lr_pilon.pdf  ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/SE02_002/final_assembly/
scp marimond@dell-node-12.mpipz.mpg.de:marimond@dell-node-12:/biodata/dep_coupland/common/Arabis_alpina_resource/annotation/annotation_methods.pdf  ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/SE02_002/final_assembly/
scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/final_assembly/annotation_methods.pdf  ~/Documents/Studium/Köln_Biological\ Science/Master\ Thesis/assemblies/SE02_002/final_assembly/
#analysis
echo "create bed"
awk '/^chr1_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr1.bed
awk '/^chr2_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr2.bed
awk '/^chr3_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr3.bed
awk '/^chr4_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr4.bed
awk '/^chr5_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr5.bed
awk '/^chr6_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr6.bed
awk '/^chr7_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr7.bed
awk '/^chr8_RagTag_pilon_pilon/{printf("%s\t0\t%s\n",$1,$2);}' SE02_002_final_assembly.fai > SE02_002_chr8.bed
echo "create individual bams"
samtools view -L SE02_002_chr1.bed -o SE02_002_chr1.bam SE02_002_pilon_long_reads.bam
samtools view -L SE02_002_chr2.bed -o SE02_002_chr2.bam SE02_002_pilon_long_reads.bam
samtools view -L SE02_002_chr3.bed -o SE02_002_chr3.bam SE02_002_pilon_long_reads.bam
samtools view -L SE02_002_chr4.bed -o SE02_002_chr4.bam SE02_002_pilon_long_reads.bam
samtools view -L SE02_002_chr5.bed -o SE02_002_chr5.bam SE02_002_pilon_long_reads.bam
samtools view -L SE02_002_chr6.bed -o SE02_002_chr6.bam SE02_002_pilon_long_reads.bam
samtools view -L SE02_002_chr7.bed -o SE02_002_chr7.bam SE02_002_pilon_long_reads.bam
samtools view -L SE02_002_chr8.bed -o SE02_002_chr8.bam SE02_002_pilon_long_reads.bam
echo "start indexing"
samtools index -b SE02_002_chr1.bam
samtools index -b SE02_002_chr2.bam
samtools index -b SE02_002_chr3.bam
samtools index -b SE02_002_chr4.bam
samtools index -b SE02_002_chr5.bam
samtools index -b SE02_002_chr6.bam
samtools index -b SE02_002_chr7.bam
samtools index -b SE02_002_chr8.bam
echo "tidy up"
mkdir bed_files
mv *bed bed_files/
echo "------------------------------------------------------------------------------------------------"
echo "create copy file"
for i in {1..8}; do
  echo "
  scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/final_assembly/SE02_002_chr${i}.bam
  scp marimond@dell-node-12.mpipz.mpg.de:/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/final_assembly/SE02_002_chr${i}.bam.bai
  " >> copy_bam_and_fasta.txt
done

echo "all done""

#fa to fasta
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' SE02_002_final_assembly.fasta | tr "\t" "\n" > SE02_002_final_assembly.fasta
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' pangenome_protein_library.fna  | tr "\t" "\n" > pangenome_protein_library.fasta
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' pangenome_transcripts_library_final.fasta  | tr "\t" "\n" > pangenome_transcripts_library_final.fa


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


cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/AT03_004/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py AT03_004_final_assembly.fasta  -o quast_merged_AT03_004 -e --fast -b &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/CZ01_005/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py CZ01_005_final_assembly.fasta  -o quast_merged_CZ01_005 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES03_014/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py ES03_014_final_assembly.fasta  -o quast_merged_ES03_014 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES18_040/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py ES18_040_final_assembly.fasta  -o quast_merged_ES18_040 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES20_028/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py ES20_028_final_assembly.fasta  -o quast_merged_ES20_028 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/ES24_001/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py ES24_001_final_assembly.fasta  -o quast_merged_ES24_001 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/NO04_015/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py NO04_015_final_assembly.fasta  -o quast_merged_NO04_015 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/PL01_004/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py PL01_004_final_assembly.fasta  -o quast_merged_PL01_004 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/PT01_005/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py PT01_005_final_assembly.fasta  -o quast_merged_PT01_005 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py SE02_002_final_assembly.fasta  -o quast_merged_SE02_002 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE06_020/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py SE06_020_final_assembly.fasta  -o quast_merged_SE06_020 -e --fast &
cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/SE07_013/final_assembly
/home/marimond/.conda/envs/bases/bin/quast.py SE07_013_final_assembly.fasta  -o quast_merged_SE07_013 -e --fast &


##########ANALYSIS FINAL##################

input_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/SC_assemblies.txt"
analysis="/netscratch/dep_coupland/grp_fulgione/male/assemblies/analysis_final_assemblies.txt"
target_character="N"
while read -r col1 col2 _; do
  cd /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/
  echo "analyse ${col1}"
  echo "${col1}" >> $analysis
  echo "MERGED CONTIGS:" >> $analysis
  echo "Busco" >> $analysis
  cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/merged_flye_
  /busco_merged_${col1}/run_brassicales_odb10/short_summary.txt >> $analysis
  echo "Quast" >> $analysis
  cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/merged_flye_hifiasm/quast_merged_${col1}/report.txt >> $analysis
  echo "PSEUDOCHROMOSOMES:" >> $analysis
  echo "Busco" >> $analysis
  cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly/busco_merged_${col1}/run_brassicales_odb10/short_summary.txt >> $analysis
  echo "Quast" >> $analysis
  cat /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/final_assembly/quast_merged_${col1}/report.txt >> $analysis
  target_character="N"
  fasta_file="/netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/pilon/pilon.fasta"
  grep -v "^>" "$fasta_file" | tr -d '\n' > sequences.txt
  awk -v char="$target_character" '{print gsub(char, "")}' sequences.txt | paste -sd+ | bc >> $analysis
  rm sequences.txt
  samtools faidx /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/pilon/pilon.fasta
  cut -f1,2 /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/pilon/pilon.fasta.fai >> $analysis
  cut -f2,3,4 /netscratch/dep_coupland/grp_fulgione/male/assemblies/${col1}/RagTag_nucmer/ragtag_output/ragtag.scaffold.confidence.txt | awk '{ total_grouping += $1; total_location += $2; total_orientation += $3; count++ } END { print "Mean Grouping Confidence:", total_grouping/count; print "Mean Location Confidence:", total_location/count; print "Mean Orientation Confidence:", total_orientation/count }' >> $analysis
  echo "--------------------------------------------------------------------------------------------" >> $analysis
done < "$input_file"




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

# Read the FASTA file and extract matching contigs
while IFS= read -r line; do
    if [[ $line =~ ^\> ]]; then
        # Check if the line is a FASTA header
        header="$line"
        contig_name=$(echo "$header" | sed 's/^>//')
        if [[ ${contigs["$contig_name"]} ]]; then
            found=1
            echo "$header"
        else
            found=0
        fi
    else
        # Output sequence if the corresponding header was found
        if [[ $found -eq 1 ]]; then
            echo "$line"
        fi
    fi
done < "$fasta_file" >> "$output_file"

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
