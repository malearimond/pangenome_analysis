#define input file
#input <- read.table("/netscratch/dep_coupland/grp_fulgione/male/assemblies/IS01_007/gene_pangenome/IS01_007_pan-gene.coverage_HEAD.coverage",header = T)
args = commandArgs(trailingOnly=TRUE)
input_path = args[1]
input <- read.table(input_path,header = T)

#define list of unique genes
unique_genes <- unique(input$GENE)

#input[which(input$GENE==unique_genes[2]),]

for(j in c(1:length(unique_genes))){
total_count <- c()
count <- 0
for(i in c(1:length(input[which(input$GENE==unique_genes[j]),][which(input[which((input$GENE==unique_genes[j])),]$DEPTH==0),]))){
 count <- count +(input[which(input$GENE==unique_genes[j]),][which(input[which((input$GENE==unique_genes[j])),]$DEPTH==0),][i,c(3)]-input[which(input$GENE==unique_genes[j]),][which(input[which((input$GENE==unique_genes[j])),]$DEPTH==0),][i,c(2)])
}
total_count <- append(total_count,(count/max(input[which(input$GENE==unique_genes[j]),3])))
}

coverage_list <- matrix(data = NA,nrow = length(unique_genes),ncol = 2)
coverage_list[,1] <- unique_genes
coverage_list[,2] <- total_count


write.table(coverage_list,file = "Coverage_genes")
