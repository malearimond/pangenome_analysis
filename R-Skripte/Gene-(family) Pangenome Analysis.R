####### ANALYSIS OF GENE-(FAMILY) PANGENOME #######
#load packages
library(ggpubr)
#install.packages("cli")
library(cli)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(VennDiagram)

#set working directory 
setwd()
###############################################################################
###### 1.) Unique gene Analysis #####
##### 1.1) Principal component analysis (PCA) of Gene PAV #########
####  1.1.1) with all samples included####
#load and show data
gene_presence_absence_matrix <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/presence_absence_matrix_0.85.tsv")
dim(gene_presence_absence_matrix)
View(gene_presence_absence_matrix)
#transform as matrix
as.matrix(gene_presence_absence_matrix)
m <- as.matrix(gene_presence_absence_matrix[, -1])    #remove first row
dim(m)
rownames(m) <- gene_presence_absence_matrix$Gene      #add gene names as rownames
n <- t(m)                                             #transpose matrix

#remove all columns which have same entry
samecol <- apply(n, 2, function(col) length(unique(col)) == 1)  
n <- n[, !samecol]
sum(is.na(n))

#compute PCA
pca_result <- prcomp(n, scale. = TRUE)
row.names(n)

#create population vector in order of rownames
pop <- c('La Palma', 'Spain', 'Poland', 'Spain', 'Branavieja', 'Swiss', 'Norway', 'Slowenia', 'Germany', 'Sweden','Norway', 'Montgarri','Sweden', 'Austria', 'Portugal','Norway', 'Spain', 'Iceland', 'Spain', 'Czech', 'Sweden' )
n <- cbind(pop, n)

# Define custom colors and shapes for the populations
colors <- c('Germany' = '#351E10', 'Norway' = '#0072B2', 'Spain' = 'darkgreen', 'Czech' = '#613613', 'Poland' = '#D55E00', 'Portugal' = 'lightgreen', 'Slowenia' = 'brown', 'Austria' = '#E69F00', 'Iceland' = '#56B4E9', 'Sweden' = 'darkblue', 'Swiss' = 'darkred','La Palma'='#DCE319FF', 'Montgarri'='#B8DE29FF', 'Branavieja'='#95D840FF')
shapes <- c('Spain' = 0, 'Norway' = 1, 'Germany' = 2, 'Sweden' = 3, 'Poland' = 4, 'Portugal' = 5, 'Iceland' = 6, 'Austria' = 7, 'Slowenia' = 8, 'Czech' = 9, 'Swiss' = 10, 'La Palma'=11, 'Montgarri'=12, 'Branavieja'=13)

#visualize PCA results
biplot <- fviz_pca_ind(
  pca_result,
  label = 'ind',  
  ggtheme = theme_minimal(),
  habillage = n[,'pop'],
  alpha.ind = 0.6,
  pointshape.ind = 1,
  title = 'PCA Plot of Gene PAV across Populations',
  xlim = NULL, , repel = TRUE
)+
  scale_color_manual(values = colors, guide = guide_legend(title = 'Populations')) +
  scale_shape_manual(values = shapes, guide = guide_legend(title = 'Populations'))+
  theme_minimal()


# Print the biplot
print(biplot)



#### 1.1.2) with low coverage samples excluded####
presence_absence_matrix_0.85 <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/presence_absence_matrix_0.85.tsv")
columns_to_remove <- c("AT03_004", "ES27_024", "CZ01_005")            #define samples to be removed
df <- presence_absence_matrix_0.85

df <- df[, !names(df) %in% columns_to_remove]                         # Remove the specified columns

# Print the updated data frame
m <- as.matrix(df[, -1])
dim(m)
rownames(m) <- presence_absence_matrix_0.85$Gene
n <- t(m)
samecol <- apply(n, 2, function(col) length(unique(col)) == 1)
n <- n[, !samecol]
sum(is.na(n))

#Compute PCA
pca_result <- prcomp(n, scale. = TRUE)

#create vector of populations and regions
pop <- c('Spain', 'Poland', 'Spain', 'Branavieja', 'Swiss', 'Norway', 'Slowenia', 'Germany', 'Sweden','Norway', 'Montgarri','Sweden', 'Portugal','Norway', 'Spain', 'Iceland', 'Spain', 'Sweden' )
Region <- c('Westeurope', 'Easteurope', 'Westeurope', 'Westeurope', 'Central-europe/Alps', 'Scandinavia', 'Central-europe/Alps','Central-europe/Alps','Scandinavia', 'Scandinavia','Westeurope', 'Scandinavia','Westeurope', 'Scandinavia', 'Westeurope', 'Scandinavia', 'Westeurope', 'Scandinavia')

#bind vectors with dataframe 
n <- cbind(pop, n)
n <- cbind(Region, n)

# Define custom colors and shapes for the 12 groups
colors <- c('Germany' = '#351E10', 'Norway' = '#0072B2', 'Spain' = 'darkgreen', 'Poland' = '#D55E00', 'Portugal' = 'lightgreen', 'Slowenia' = 'brown', 'Iceland' = '#56B4E9', 'Sweden' = 'darkblue', 'Swiss' = 'darkred', 'Montgarri'='#B8DE29FF', 'Branavieja'='#95D840FF')
shapes <- c('Spain' = 0, 'Norway' = 1, 'Germany' = 2, 'Sweden' = 3, 'Poland' = 4, 'Portugal' = 5, 'Iceland' = 6, 'Slowenia' = 8,  'Swiss' = 10,  'Montgarri'=12, 'Branavieja'=13)

#visualize PCA
biplot <- fviz_pca_ind(
  pca_result,
  label = 'ind',  
  ggtheme = theme_minimal(),
  habillage = n[,'pop'],
  alpha.ind = 0.6,
  pointshape.ind = 1,
  title = 'PCA Plot of Gene PAV across Populations',
  xlim = NULL, , repel = TRUE
)+
  scale_color_manual(values = colors, guide = guide_legend(title = 'Populations')) +
  scale_shape_manual(values = shapes, guide = guide_legend(title = 'Populations'))+
  theme_minimal()


# Print the biplot
print(biplot)


summary(pca_result)


#### 1.2) Create Venn Diagram of Gene PAV####
presence_absence_matrix_0.85 <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/presence_absence_matrix_0.85.tsv", row.names = 1)
columns_to_remove <- c("AT03_004", "ES27_024", "CZ01_005")
df <- presence_absence_matrix_0.85
# Remove the specified columns
df <- df[, !names(df) %in% columns_to_remove]
colnames(df)

## Subsample to specific regions
#Scandinavian 
Scandinavia <- df[, c("NO01_005", "NO04_015", "SE02_002", "NO02_001", "IS01_007", "SE06_020","SE07_013")]

Sums_of_Scandinavia<-rowSums(Scandinavia)                                   #Sum up all Gene present counts
Scandinavia$Scandinavia_Sum <- Sums_of_Scandinavia                          #Add vector to scandinavian data_frame
genes_scandinavia <- row.names(Scandinavia[Scandinavia$Scandinavia_Sum==7,])#Filter genes which are present in all scandinavian ecotypes
write.table(genes_scandinavia, "Scandinavia_genes_only.txt")                #Save all gene names present in scandinavia

#Westeurope (Spain/Portugal)
Spain <- df[, c("ES03_014","ES04_014", "ES18_040","ES20_028", "ES06_006","ES24_001","PT01_005")]
Sums_of_Spain<-rowSums(Spain)
Spain$Spain_Sum <- Sums_of_Spain
genes_spain <- row.names(Spain[Spain$Spain_Sum==7,])
write.table(genes_spain, "Spain_genes_only.txt")

#Central_europe
Central_eupope <- df[, c("SI01_005", "DE01_001","CH05_002")]
Sums_of_Central_eupope<-rowSums(Central_eupope)
Central_eupope$Central_eupope_Sum <- Sums_of_Central_eupope
genes_central_europe <- row.names(Central_eupope[Central_eupope$Central_eupope_Sum==3,])
write.table(genes_central_europe, "Central_europe_genes_only.txt")

#East_europe (Poland)
genes_poland <- row.names(df[df$PL01_004==1,])
write.table(genes_poland, "Poland_genes_only.txt")

# Create a list of gene vectors 

gene_list <- list(
  vector1 = genes_scandinavia,
  vector2 = genes_central_europe,
  vector3 = genes_spain,
  vector4 = genes_poland
)


# Initialize a list to store unique non-overlapping genes
unique_genes_list <- list()

# Loop through each vector of genes
for (i in 1:length(gene_list)) {
  current_vector <- gene_list[[i]]
  non_overlapping_genes <- current_vector
  
  # Check for overlap with other vectors and remove overlapping genes
  for (j in 1:length(gene_list)) {
    if (i != j) {
      non_overlapping_genes <- setdiff(non_overlapping_genes, gene_list[[j]])
    }
  }
  
  unique_genes_list[[i]] <- non_overlapping_genes
}


#Unique genes per region
genes_scandinavia_unique <- unique_genes_list[[1]]
genes_central_europe_unique <- unique_genes_list[[2]]
genes_spain_unique <- unique_genes_list[[3]]
genes_poland_unique <- unique_genes_list[[4]]

#write table for unqiue genes
write.table(genes_scandinavia_unique, file="scandinavia_genes_unique.txt", row.names = F, col.names = F)
write.table(genes_spain_unique, file="spain_genes_unique.txt", row.names = F, col.names = F)
write.table(genes_poland_unique, file="poland_genes_unique.txt", row.names = F, col.names = F)
write.table(genes_central_europe_unique, file="central_europe_genes_unique.txt", row.names = F, col.names = F)


#VISUALIZE VENN DIAGRAM
#Define function for venn_diagramm
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

#venn.diagram(x, filename = "venn-4-dimensions.png")

cbPalette_9 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")
venn_object <- venn.diagram(gene_list,
                            filename = "Venn_Arabis_alpina.png",
                            category.names = c("Poland" , "Spain/Portugal" , "Scandinavia", "Central-Europe"),
                            # Circles
                            lwd = 2,
                            lty = 'blank',
                            fill = c("#999999", "#E69F00", "#56B4E9", "#D55E00"),
                            # Numbers
                            cex = .9,
                            fontface = "italic",
                            # Set names
                            cat.cex = 1,
                            cat.fontface = "bold",
                            cat.default.pos = "outer",
                            cat.dist = c(0.055, 0.055, 0.1, 0.1),
                            main = "Pangenome of Arabis alpina",
                            main.cex=1.5
                            )

display_venn(
  gene_list,
  category.names = c("Poland" , "Spain/Portugal" , "Scandinavia", "Central-Europe"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#999999", "#E69F00", "#56B4E9", "#D55E00"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1),
  main = "Pangenome of Arabis alpina",
  main.cex=1.5
)

#### 1.3) Compute Gene Pangenome Categories ####
#Core and dispensable genes/cloud genes
presence_absence_matrix_0.85 <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/presence_absence_matrix_0.85.tsv", row.names = 1)
columns_to_remove <- c("AT03_004", "ES27_024", "CZ01_005")
df <- presence_absence_matrix_0.85
# Remove the specified columns
df <- df[, !names(df) %in% columns_to_remove]
Sums_of_present_genes<-rowSums(df)
df$Sum <- Sums_of_present_genes

#define categories
Core <- length(df[df$Sum==18,]$Sum)
private <- length(df[df$Sum==1,]$Sum)
soft_core <- length(df[df$Sum > 15 & df$Sum < 18,]$Sum)
dispensable <- length(df[df$Sum > 1 & df$Sum < 16 ,]$Sum)
sum(Core, soft_core, dispensable, private)                          #valid all genes are captured?

#extract gene names per category
core_genes_list <-row.names(df[df$Sum==18,])
private_genes_list <- row.names(df[df$Sum==1,])
soft_core_genes_list <- row.names(df[df$Sum > 15 & df$Sum < 18,])
dispensable_genes_list <- row.names(df[df$Sum > 1 & df$Sum < 16 ,])

#save gene names per category
write.table(core_genes_list, "core_genes.txt", row.names = F, col.names = F)
write.table(private_genes_list, "private_genes.txt", row.names = F, col.names = F)
write.table(soft_core_genes_list, "soft_core_genes.txt", row.names = F, col.names = F)
write.table(dispensable_genes_list, "dispensable_genes.txt", row.names = F, col.names = F)

#Add column with corresponding pangenome category
df$Class[df$Sum == 18] <- "core"
df$Class[df$Sum == 1] <- "private"
df$Class[df$Sum > 15 & df$Sum < 18] <- "soft_core"
df$Class[df$Sum > 1 & df$Sum < 16] <- "dispensable"

#Visualize Pangenome properties
hist <-ggplot(data = subset(df, !is.na(Class)), aes(x = Sum, fill = Class, na.rm  =TRUE)) +
  geom_histogram(bins = 18, binwidth = 1, color = "black")+
  scale_fill_manual(values=cbPalette_9) + theme_minimal()+
  labs(title = paste("Gene Variation Pangenome unique genes") , x="Number of accession", y = "Frequency")+
  ylim(0, 18000)

#Create circular plot (first write table with relative frequencies)
Realtive_frequencies <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/Realtive_frequencies.txt")
Realtive_frequencies$Relative.Frequency
Realtive_frequencies$Prozent <- paste(round(Realtive_frequencies$Frequency/sum(Realtive_frequencies$Frequency)*100,2), "%")

circle <- ggplot(Realtive_frequencies, aes(x = "", y=Frequency, fill=Pangenome.category))+
  geom_bar(stat="identity", color="black")+
  coord_polar("y")+
  theme_void()+
  scale_fill_manual(values=cbPalette_9)+
  geom_label(aes(label=Prozent),
             position = position_stack(vjust = 0.5),
             color=c("black", "black", "black", "black"),
             label.size = 0,
             size=3, 
             show.legend = FALSE)+
  theme(legend.position = "none")

#Arrange both plots
ggarrange(hist, circle)



#########################################################################################
######## 2.) Orthogroup Analysis #############

##### 2.1) Principal component analysis of Gene-family PAV #####
Orthogroups <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/Orthogroups.GeneCount.tsv")

df <- Orthogroups
columns_to_remove <- c("AT03_004_proteins", "ES27_024_proteins", "CZ01_005_proteins", "Total",
                       "Alyrata_384_v2.1.protein","Aalpina_v5.1_proteins", "Athaliana_167_protein", "Crubella_474_v1.1.protein")

# Remove the specified columns
df <- df[, !names(df) %in% columns_to_remove]
colnames(df)
colnames(df)  <- c("Orthogroup","CH05_002", "DE01_001", "ES03_014", "ES04_014", "ES06_006","ES18_040", "ES20_028", "ES24_001", "IS01_007", "NO01_005","NO02_001", "NO04_015", "PL01_004", "PT01_005", "SE02_002","SE06_020", "SE07_013", "SI01_005")

m <- as.matrix(df[, -1])

rownames(m) <- Orthogroups$Gene
head(m)
n <- t(m)
samecol <- apply(n, 2, function(col) length(unique(col)) == 1)
n <- n[, !samecol]
sum(is.na(n))

#Compute PCA
pca_result <- prcomp(n, scale. = TRUE)
row.names(n)

pop <- c('Swiss','Germany','Spain',  'Spain', 'Spain',  'Branavieja', 'Montgarri',  'Spain','Iceland', 'Norway','Norway', 'Norway','Poland','Portugal', 'Sweden', 'Sweden', 'Sweden', 'Slowenia' )
Region <- c('Westeurope', 'Easteurope', 'Westeurope', 'Westeurope', 'Central-europe/Alps', 'Scandinavia', 'Central-europe/Alps','Central-europe/Alps','Scandinavia', 'Scandinavia','Westeurope', 'Scandinavia','Westeurope', 'Scandinavia', 'Westeurope', 'Scandinavia', 'Westeurope', 'Scandinavia')
length(pop)
length(Region)
n <- cbind(pop, n)
n <- cbind(Region, n)


# Define custom colors and shapes for the orthogroups
cbPalette_9 <- c("darkgrey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")

colors <- c('Germany' = '#351E10', 'Norway' = '#0072B2', 'Spain' = '#D55E00', 'Poland' ="#999999" , 'Portugal' = 'orange', 'Slowenia' = 'brown', 'Iceland' = '#56B4E9', 'Sweden' = 'darkblue', 'Swiss' = 'darkred', 'Montgarri'='#E69F00', 'Branavieja'= "#F0E442")
shapes <- c('Spain' = 0, 'Norway' = 1, 'Germany' = 2, 'Sweden' = 3, 'Poland' = 4, 'Portugal' = 5, 'Iceland' = 6, 'Slowenia' = 8,  'Swiss' = 10,  'Montgarri'=12, 'Branavieja'=13)
biplot <- fviz_pca_ind(
  pca_result,
  label = 'ind',  
  ggtheme = theme_minimal(),
  habillage = n[,'pop'],
  alpha.ind = 0.6,
  pointshape.ind = 1,
  title = 'PCA Plot of Gene family PAV across Populations',
  ,  repel = TRUE
)+
  scale_color_manual(values = colors, guide = guide_legend(title = 'Populations')) +
  scale_shape_manual(values = shapes, guide = guide_legend(title = 'Populations'))+
  theme_minimal()


# Print the biplot
print(biplot)
summary(pca_result)


#### 2.2) Compute Venn Diagram for Orthogroups ####

#read table
Orthogroups <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/Orthogroups.GeneCount.tsv", row.names = 1)

#Remove columns 
columns_to_remove <- c("AT03_004_proteins", "ES27_024_proteins", "CZ01_005_proteins", "Total",
                       "Alyrata_384_v2.1.protein","Aalpina_v5.1_proteins", "Athaliana_167_protein", "Crubella_474_v1.1.protein")
columns_to_remove <- c("Total","CZ01_005_proteins",
                       "Alyrata_384_v2.1.protein","Aalpina_v5.1_proteins", "Athaliana_167_protein", "Crubella_474_v1.1.protein")
df <- Orthogroups
df[df > 1] <- 1
# Remove the specified columns
df <- df[, !names(df) %in% columns_to_remove]

#Create Groups/Regions of Venn Diagram
#Scandinavian 
Scandinavia <- df[, c("NO01_005_proteins", "NO04_015_proteins", "SE02_002_proteins", "NO02_001_proteins", "IS01_007_proteins", "SE06_020_proteins","SE07_013_proteins")]
Sums_of_Scandinavia<-rowSums(Scandinavia)
Scandinavia$Scandinavia_Sum <- Sums_of_Scandinavia
genes_scandinavia <- row.names(Scandinavia[Scandinavia$Scandinavia_Sum==7,])

#create Westerneurope
Spain <- df[, c("ES27_024_proteins","ES03_014_proteins","ES04_014_proteins", "ES18_040_proteins","ES20_028_proteins", "ES06_006_proteins","ES24_001_proteins","PT01_005_proteins")]
Sums_of_Spain<-rowSums(Spain)
Spain$Spain_Sum <- Sums_of_Spain
genes_spain <- row.names(Spain[Spain$Spain_Sum==8,])

#central_europe
Central_eupope <- df[, c("AT03_004_proteins",  "SI01_005_proteins", "DE01_001_proteins","CH05_002_proteins")]
Sums_of_Central_eupope<-rowSums(Central_eupope)
Central_eupope$Central_eupope_Sum <- Sums_of_Central_eupope
genes_central_europe <- row.names(Central_eupope[Central_eupope$Central_eupope_Sum==4,])

#Easteurope
genes_poland <- row.names(df[df$PL01_004_proteins==1,])
east_eupope <- subset(df, select = PL01_004_proteins)
Poland <- row.names(east_eupope[east_eupope$PL01_004>0,])

#Summarize all orthogroup sets into one list
x <- list(A = genes_poland,B =  genes_spain,C = genes_scandinavia, D =genes_central_europe)



#Define function for venn_diagramm
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

cbPalette_9 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")
venn_object <- venn.diagram(x,
                            filename = "Venn_Arabis_alpina_Orthogroups.png",
                            category.names = c("Poland" , "Spain/Portugal" , "Scandinavia", "Central-Europe"),
                            # Circles
                            lwd = 2,
                            lty = 'blank',
                            fill = c("#999999", "#E69F00", "#56B4E9", "#D55E00"),
                            # Numbers
                            cex = .9,
                            fontface = "italic",
                            # Set names
                            cat.cex = 1,
                            cat.fontface = "bold",
                            cat.default.pos = "outer",
                            cat.dist = c(0.055, 0.055, 0.1, 0.1),
                            main = "Pangenome of Arabis alpina Orthogroups",
                            main.cex=1.5
)


display_venn(
  x,
  category.names = c("Poland" , "Spain/Portugal" , "Scandinavia", "Central-Europe"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#999999", "#E69F00", "#56B4E9", "#D55E00"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1),
  main = "Pangenome of Arabis alpina Orthogroups",
  main.cex=1.5
)




##### 2.3) Orthogroup families (core/dispensable), circle#####
Orthogroups <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/Orthogroups.GeneCount.tsv", row.names = 1)

columns_to_remove <- c("AT03_004_proteins", "ES27_024_proteins", "CZ01_005_proteins",
                       "Alyrata_384_v2.1.protein","Aalpina_v5.1_proteins", "Athaliana_167_protein", "Crubella_474_v1.1.protein", "Total")
df <- Orthogroups
# Remove the specified columns
df <- df[, !names(df) %in% columns_to_remove]

# Set values greater than 1 to 1
df[df > 1] <- 1

# Print the updated data frame
m <- as.matrix(df[, -1])

#set all values of orthogroups to 0 and 1 to create matrix
m[m > 1] <- 1
m
write.table(m, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/0_1_matrix_gene_families.txt", sep = "", quote = FALSE, row.names=F, col.names=F)


#Determine cor/private Gene families
Sums_of_present_genes<-rowSums(df)
df$Sum <- Sums_of_present_genes
# Create a new dataframe excluding rows where df$Sum == 0
df_filtered_gene_families <- df[df$Sum != 0, ]
Core <- length(df_filtered_gene_families[df_filtered_gene_families$Sum==18,]$Sum)
private <- length(df_filtered_gene_families[df_filtered_gene_families$Sum==1,]$Sum)
soft_core <- length(df_filtered_gene_families[df_filtered_gene_families$Sum > 15 & df_filtered_gene_families$Sum < 18,]$Sum)
dispensable <- length(df_filtered_gene_families[df_filtered_gene_families$Sum > 1 & df_filtered_gene_families$Sum < 16 ,]$Sum)
sum(Core, soft_core, dispensable, private)

#Create Plots for Gene families#
df_filtered_gene_families$Class[df_filtered_gene_families$Sum == 18] <- "core"
df_filtered_gene_families$Class[df_filtered_gene_families$Sum==1] <- "private"
df_filtered_gene_families$Class[df_filtered_gene_families$Sum > 15 & df_filtered_gene_families$Sum < 18] <- "soft_core"
df_filtered_gene_families$Class[df_filtered_gene_families$Sum > 1 & df_filtered_gene_families$Sum < 16] <- "dispensable"


cbPalette_9 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


hist_gene_families <-ggplot(data = subset(df_filtered_gene_families, !is.na(Class)), aes(x = Sum, fill = Class, na.rm  =TRUE)) +
  geom_histogram(bins = 18, binwidth = 1, color = "black")+
  scale_fill_manual(values=cbPalette_9) + theme_minimal()+
  labs(title = paste("Gene Variation Pangenome Gene families") , x="Number of accession", y = "Frequency")+
  ylim(0, 18000)

#Circulat plot
Realtive_frequencies <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/Realtive_frequencies_gene_families.txt")
Realtive_frequencies$Relative.Frequency
Realtive_frequencies$Prozent <- paste(round(Realtive_frequencies$Frequency/sum(Realtive_frequencies$Frequency)*100,2), "%")

circle_gene_families <- ggplot(Realtive_frequencies, aes(x = "", y=Frequency, fill=Pangenome.category))+
  geom_bar(stat="identity", color="black")+
  coord_polar("y")+
  theme_void()+
  scale_fill_manual(values=cbPalette_9)+
  geom_label(aes(label=Prozent),
             position = position_stack(vjust = 0.5),
             color=c("black", "black", "black", "black"),
             label.size = 0,
             size=3, 
             show.legend = FALSE)+
  theme(legend.position = "none")


ggarrange(hist_gene_families, circle_gene_families)


################## 3.) COMPARE Annotation methods and coverage correlation############
Gene_presence_comparison <- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/Gene_presence_comparison.txt")

Gene_presence_comparison$Mapping
Gene_presence_comparison$Liftover

ggplot(data = Gene_presence_comparison, aes(x=Coverage,y=Liftover, colour=Method ))+
  geom_point(size=2)+
  scale_colour_manual(values=cbPalette_9)+
  ylab("Number of genes present")+
  ggtitle("Observed number of genes depending on \"annotation\" method")+
  theme_minimal()
