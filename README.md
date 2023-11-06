# MASTER THESIS
Github repository consists of all scripts which were used to perform de novo assembly for 21 Hifi raw read sets. Additionally genome annotation and gene-/family pangenome were constructed and further analyzed in RStudio.

# Assembly workflow
Whole assembly workflow can be found in Assembly.sh. Two assembly tools (hifiasm and Flye) were used to perform the assembly. Afterwards they were merged using quickmerge and depending on coverage purged using purge_dups. Scaffolding and Gapclosing was performed with RagTag and LR_Gapcloser. Evaluation was conducted using BUSCO and QUAST.

# Construction of a Gene-family Pangenome
Gene-family pangenome was constructed using revealed genes via Map-to-Pan strategy which were clustered with OrthoFinder. As outgroups served A. thaliana, A.lyrata, C. rubella (and additionally A. alpina v5.1).


