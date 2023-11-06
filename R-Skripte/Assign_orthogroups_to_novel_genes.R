#Get Ortholog name for each contig
library(dplyr)

orthologs <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/Orthogroups.tsv",sep ="\t",header = T)

subset_orthologs <- subset_orthologs %>%
  mutate(across(everything(), ~na_if(., "")))
View(subset_orthologs)

##########################################
# Create an empty data frame to store the results
result_df <- data.frame(Orthogroup = character(0), Contig = character(0), Description= character(0), GeneID = character(0), Ortholog_species=character(0), stringsAsFactors = FALSE)
number_of_novel_genes <- 0
#for (i in 1:nrow(subset_orthologs))
###Test 1:2
for (i in 1:nrow(subset_orthologs)) {
# Get unique contig names from the first row
row_entries <- as.character(subset_orthologs[i, ])
entries <- unlist(strsplit(row_entries, " "))
contig_names <- unique(grep("^_contig_\\S+\\.p[123]", entries, value = TRUE))
protein_col <- NULL
Ortholog_species <- NULL
# Iterate through each row
  
  contig <- contig_names
  gene_id <- NA
  Orthogroup <- rep(subset_orthologs[i, "Orthogroup"], length(contig))
  
  # Check the protein columns in the order you specified
  if (!is.na(subset_orthologs[i, "Aalpina_v5.1_proteins"])) {
    protein_col <- rep(subset_orthologs[i, "Aalpina_v5.1_proteins"], length(contig))
    entries_pro <- unlist(strsplit(protein_col, " "))
    gene_id <- rep(entries_pro[1], length(contig))
    ortholog_species <- rep("Arabis alpina", length(contig))
  } else if (!is.na(subset_orthologs[i, "Athaliana_167_protein"])) {
    protein_col <- rep(subset_orthologs[i, "Athaliana_167_protein"], length(contig))
    entries_pro <- unlist(strsplit(protein_col, " "))
    gene_id <- rep(entries_pro[1], length(contig))
    ortholog_species <- rep("Arabidopsis thaliana", length(contig))
  } else if (!is.na(subset_orthologs[i, "Alyrata_384_v2.1.protein"])) {
    protein_col <- rep(subset_orthologs[i, "Alyrata_384_v2.1.protein"], length(contig))
    entries_pro <- unlist(strsplit(protein_col, " "))
    gene_id <- rep(entries_pro[1], length(contig))
    ortholog_species <- rep("Arabidopsis lyrata", length(contig))
  } else if (!is.na(subset_orthologs[i, "Crubella_474_v1.1.protein"])) {
    protein_col <- rep(subset_orthologs[i, "Crubella_474_v1.1.protein"], length(contig))
    entries_pro <- unlist(strsplit(protein_col, " "))
    gene_id <- rep(entries_pro[1], length(contig))
    ortholog_species <- rep("Crubella", length(contig))
  } else {
    protein_col <- rep("NA", length(contig))
    gene_id <- protein_col
    ortholog_species <- rep("NA", length(contig))
  }
  
  # Add the information to the result data frame
  number_of_novel_genes <- length(contig) + sum(number_of_novel_genes)
  result_df <- rbind(result_df, data.frame(Orthogroup=Orthogroup, Contig = contig, Description= protein_col,GeneID = gene_id, Ortholog_species = ortholog_species))
}

View(result_df)

setwd("~/Studium/Köln_Biological Science/Master Thesis/Results/")
# Save the result data frame to a file if needed
write.table(result_df, "Orthogroups_of_novel_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#total:
length(result_df$Orthogroup) #6787
length(result_df$Orthogroup[result_df$Ortholog_species=="Arabis alpina"]) #2381
length(result_df$Orthogroup[result_df$Ortholog_species=="NA"]) #2984
length(result_df$Orthogroup[result_df$Ortholog_species=="Arabidopsis thaliana"]) #987
length(result_df$Orthogroup[result_df$Ortholog_species=="Arabidopsis lyrata"]) #271
length(result_df$Orthogroup[result_df$Ortholog_species=="Crubella"]) #164


results_only_alpina <- result_df[result_df$Ortholog_species=="Arabis alpina",]
write.table(results_only_alpina, "Orthogroups_of_novel_genes_only_alpina.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#########
core_genes <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/core_genes.txt")
softcore_genes <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/soft_core_genes.txt")
disp_genes <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/dispensable_genes.txt")
private_genes <- read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/private_genes.txt")
head(core_genes)
length(which(core_genes$V1 %in% results_only_alpina$Contig))
length(which(softcore_genes$V1 %in% results_only_alpina$Contig))
length(which(disp_genes$V1 %in% results_only_alpina$Contig))
length(which(private_genes$V1 %in% results_only_alpina$Contig))


###CORE GENES###
results_only_alpina$pangenome_category <- rep("NA", length(row.names(results_only_alpina)))
old_contig_name <- NULL
core_genes_new <- core_genes
core_position <- which(core_genes$V1 %in% results_only_alpina$Contig)
for (i in core_position) {
  old_contig_name <- core_genes$V1[i]
  core_genes_new$V1[i] <- results_only_alpina$GeneID[results_only_alpina$Contig==old_contig_name]
  
  results_only_alpina$pangenome_category[results_only_alpina$Contig==old_contig_name] <- as.character("core")
  
  print(i)
}

length(which(core_genes_new$V1 %in% results_only_alpina$Contig))


### SOFT CORE ###
old_contig_name <- NULL
soft_core_genes_new <- softcore_genes
position <- which(softcore_genes$V1 %in% results_only_alpina$Contig)
for (i in position) {
  old_contig_name <- softcore_genes$V1[i]
  soft_core_genes_new$V1[i] <- results_only_alpina$GeneID[results_only_alpina$Contig==old_contig_name]
  
  results_only_alpina$pangenome_category[results_only_alpina$Contig==old_contig_name] <- as.character("softcore")
  
  print(i)
}

length(which(soft_core_genes_new$V1 %in% results_only_alpina$Contig))

### DISPENSABLE ###
old_contig_name <- NULL
dispensable_genes_new <- disp_genes
position <- which(disp_genes$V1 %in% results_only_alpina$Contig)
for (i in position) {
  old_contig_name <- disp_genes$V1[i]
  dispensable_genes_new$V1[i] <- results_only_alpina$GeneID[results_only_alpina$Contig==old_contig_name]
  
  #results_only_alpina$pangenome_category[results_only_alpina$Contig==old_contig_name] <- as.character("softcore")
  
  print(i)
}

length(which(dispensable_genes_new$V1 %in% results_only_alpina$Contig))

### PRIVATE ###
old_contig_name <- NULL
private_genes_new <- private_genes
position <- which(private_genes$V1 %in% results_only_alpina$Contig)
for (i in position) {
  old_contig_name <- private_genes$V1[i]
  private_genes_new$V1[i] <- results_only_alpina$GeneID[results_only_alpina$Contig==old_contig_name]
  
  results_only_alpina$pangenome_category[results_only_alpina$Contig==old_contig_name] <- as.character("softcore")
  
  print(i)
}

length(which(private_genes_new$V1 %in% results_only_alpina$Contig))


####################
# Assuming you have the dplyr package installed, load it
library(dplyr)

# Keep only the unique values in the V1 column
core_genes_new <- core_genes_new %>% distinct(V1)
soft_core_genes_new <- soft_core_genes_new %>% distinct(V1)
dispensable_genes_new <- dispensable_genes_new %>% distinct(V1)
private_genes_new <- private_genes_new %>% distinct(V1)
# Remove genes in core_genes_new that are present in soft_core_genes_new


### make unique gene sets ###
core_soft_core_overlap <- c()
for (f in 1:length(which(core_genes_new$V1%in%soft_core_genes_new$V1))) {
  for (e in which(core_genes_new$V1%in%soft_core_genes_new$V1)) {
  print(core_genes_new$V1[e])
  core_soft_core_overlap[f] <- core_genes_new$V1[e]
  }
}
core_genes_new <- anti_join(core_genes_new, soft_core_genes_new, by = "V1")
length(core_genes_new$V1) #13489

#dispens
core_dispens_overlap <- c()
for (f in 1:length(which(core_genes_new$V1%in% dispensable_genes_new$V1))) {
  for (e in which(core_genes_new$V1 %in% dispensable_genes_new$V1)) {
    print(core_genes_new$V1[e])
    core_dispens_overlap[f] <- core_genes_new$V1[e]
  }
}
core_genes_new <- anti_join(core_genes_new, dispensable_genes_new, by = "V1")
length(core_genes_new$V1) #13419

core_genes_new$V1[598]

#private
core_private_overlap <- c()
for (f in 1:length(which(core_genes_new$V1%in% private_genes_new$V1))) {
  for (e in which(core_genes_new$V1 %in% private_genes_new$V1)) {
    print(core_genes_new$V1[e])
    core_private_overlap[f] <- core_genes_new$V1[e]
  }
}
core_genes_new <- anti_join(core_genes_new, private_genes_new, by = "V1")
length(core_genes_new$V1) #13379


#soft_disp
soft_disp_overlap <- c()
for (f in 1:length(which(soft_core_genes_new$V1%in%dispensable_genes_new$V1))) {
  for (e in which(soft_core_genes_new$V1%in%dispensable_genes_new$V1)) {
    print(soft_core_genes_new$V1[e])
    soft_disp_overlap[f] <- soft_core_genes_new$V1[e]
  }
}
soft_core_genes_new <- anti_join(soft_core_genes_new, dispensable_genes_new, by = "V1")
length(soft_core_genes_new$V1) #3832

#soft_priv
soft_private_overlap <- c()
for (f in 1:length(which(soft_core_genes_new$V1%in%private_genes_new$V1))) {
  for (e in which(soft_core_genes_new$V1%in%private_genes_new$V1)) {
    print(soft_core_genes_new$V1[e])
    soft_private_overlap[f] <- soft_core_genes_new$V1[e]
  }
}
soft_core_genes_new <- anti_join(soft_core_genes_new, private_genes_new, by = "V1")
length(soft_core_genes_new$V1) #3819

#disp_priv
disp_private_overlap <- c()
for (f in 1:length(which(dispensable_genes_new$V1%in%private_genes_new$V1))) {
  for (e in which(dispensable_genes_new$V1%in%private_genes_new$V1)) {
    print(dispensable_genes_new$V1[e])
    disp_private_overlap[f] <- dispensable_genes_new$V1[e]
  }
}
length(dispensable_genes_new$V1) #8758

private_genes_new <- anti_join(private_genes_new,dispensable_genes_new, by = "V1")
length(core_genes_new$V1)/(length(private_genes_new$V1)+length(core_genes_new$V1)+length(soft_core_genes_new$V1)+length(dispensable_genes_new$V1))*100

getwd()
write.table(soft_core_genes_new, "soft_core_genes_annotated.txt", row.names = F, col.names = F)
write.table(core_genes_new, "core_genes_annotated.txt", row.names = F, col.names = F)
write.table(dispensable_genes_new, "dispensable_genes_annotated.txt", row.names = F, col.names = F)
write.table(private_genes_new, "private_genes_annotated.txt", row.names = F, col.names = F)


########
#GENES
#MAF FAMILY GENES / PEP GENES

