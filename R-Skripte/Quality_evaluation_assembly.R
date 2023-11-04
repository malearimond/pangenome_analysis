##### Quality evaluation of different Assemblies and post-processing tools ####

#load packages
library(ggplot2)
library(ggpubr)


#set working directory 
directory = "~/Studium/Köln_Biological Science/Master Thesis/assemblies/"
setwd(directory)

#read in data
Quality_overview<- read.delim("~/Studium/Köln_Biological Science/Master Thesis/Results/Quality_overview_update.txt")

#### 1) Normalize and scale values between 1-10####
categories <- unique(Quality_overview$Measure)

# Create a new data frame to store normalized values
normalized_data <- data.frame()

# Loop through each category
for (category in categories) {
  subset_data <- Quality_overview[Quality_overview$Measure == category, ]  # Subset data for the current category
  
  # Calculate the minimum and maximum absolute values for the current category
  min_value <- min(subset_data$Raw_value)
  max_value <- max(subset_data$Raw_value)
  
  scale_min <- 1
  scale_max <- 10
  
  # Adjust the scale_min and scale_max for specific categories
  if (category %in% c("Complete and duplicated BUSCOs", "Number of contigs", "CPU")) {
    scale_min <- 10
    scale_max <- 1
  }
  
  normalized_values <- scale_min + (scale_max - scale_min) * ((subset_data$Raw_value - min_value) / (max_value - min_value))
  
  # Add the normalized values to the new data frame
  subset_data$NormalizedValue <- normalized_values
  normalized_data <- rbind(normalized_data, subset_data)
}

ordered_data <- normalized_data[order(normalized_data$Sample), ]
ordered_data <- ordered_data[order(ordered_data$ID), ]
View(ordered_data)

# Remove the unnecessary columns
subsetted_data <- subset(ordered_data, select = -Adjusted_Score)
subsetted_data <- subset(subsetted_data, select = -X)

#save normalized data frame
write.csv(subsetted_data, "evaluation.csv")
write.table(subsetted_data, "evaluation.txt")


#rename again into original name
Quality_overview <- ordered_data


############ 2.) NORMALIZE SC ##############
SC_Quality<- Quality_overview[Quality_overview$Mating=="SC",]

SC_categories <- unique(SC_Quality$Measure)

# Create a new data frame to store normalized values
normalized_data_SC <- data.frame()

# Loop through each category
for (category in SC_categories) {
  subset_data <- SC_Quality[SC_Quality$Measure == category, ]  # Subset data for the current category
  
  # Calculate the minimum and maximum absolute values for the current category
  min_value <- min(subset_data$Raw_value)
  max_value <- max(subset_data$Raw_value)
  
  scale_min <- 1
  scale_max <- 10
  
  # Adjust the scale_min and scale_max for specific categories
  if (category %in% c("Complete and duplicated BUSCOs", "Number of contigs", "CPU")) {
    scale_min <- 10
    scale_max <- 1
  }
  
  normalized_values <- scale_min + (scale_max - scale_min) * ((subset_data$Raw_value - min_value) / (max_value - min_value))
  
  # Add the normalized values to the new data frame
  subset_data$NormalizedValue <- normalized_values
  normalized_data_SC <- rbind(normalized_data_SC, subset_data)
}

ordered_data_SC <- normalized_data_SC[order(normalized_data_SC$Sample), ]
ordered_data_SC <- ordered_data_SC[order(ordered_data_SC$ID), ]
SC_Quality <- ordered_data_SC
SC_Quality <- SC_Quality[SC_Quality$Measure != "CPU", ]

############# NORMALIZE SI #############
SI_Quality<- Quality_overview[Quality_overview$Mating=="SI",]

SI_categories <- unique(SI_Quality$Measure)

# Create a new data frame to store normalized values
normalized_data_SI <- data.frame()

# Loop through each category
for (category in SI_categories) {
  subset_data <- SI_Quality[SI_Quality$Measure == category, ]  # Subset data for the current category
  
  # Calculate the minimum and maximum absolute values for the current category
  min_value <- min(subset_data$Raw_value)
  max_value <- max(subset_data$Raw_value)
  
  scale_min <- 1
  scale_max <- 10
  
  # Adjust the scale_min and scale_max for specific categories
  if (category %in% c("Complete and duplicated BUSCOs", "Number of contigs", "CPU")) {
    scale_min <- 10
    scale_max <- 1
  }
  
  normalized_values <- scale_min + (scale_max - scale_min) * ((subset_data$Raw_value - min_value) / (max_value - min_value))
  
  # Add the normalized values to the new data frame
  subset_data$NormalizedValue <- normalized_values
  normalized_data_SI <- rbind(normalized_data_SI, subset_data)
}

ordered_data_SI <- normalized_data_SI[order(normalized_data_SI$Sample), ]
ordered_data_SI <- ordered_data_SI[order(ordered_data_SI$ID), ]
SI_Quality <- ordered_data_SI
SI_Quality <- SI_Quality[SI_Quality$Measure != "CPU", ]

###################### 3) Create BOXPLOTS (WITHOUT CPU) for SI and SC sets #######

View(Quality_overview)
cbPalette_9 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","white")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


SI <- ggplot(SI_Quality, aes(x=NormalizedValue, y=reorder(Sample, NormalizedValue), fill=Sample)) +
  geom_boxplot()+
  scale_fill_manual(values=cbPalette_9) + theme_minimal()+
  labs(title = paste("SI: Assembly evaluation scores") , x="Quality Score", y = "Assembly Tool/Step")+
  theme(legend.position = "none")
SC <- ggplot(SC_Quality, aes(x=NormalizedValue, y=reorder(Sample, NormalizedValue), fill=Sample)) +
  geom_boxplot()+
  scale_fill_manual(values=cbPalette) + theme_minimal()+
  labs(title="SC: Assembly evaluation scores",x="Quality Score", y = "Assembly Tool/Step")+
  theme(legend.position = "none")


##### 4) Create Circular Plots #####
SI_Quality<- Quality_overview[Quality_overview$Mating=="SI",]
SC_Quality<- Quality_overview[Quality_overview$Mating=="SC",]

#Split for unique sets analyzed 
IT_Quality <- SI_Quality[SI_Quality$ID=="IT20_007",]
NO_Quality <- SC_Quality[SC_Quality$ID=="NO02_001",]
SE_Quality <- SC_Quality[SC_Quality$ID=="SE02_002",]
GR_Quality <- SI_Quality[SI_Quality$ID=="GR01_002",]

SC_polygon <- ggplot(NO_Quality, aes(x = as.factor(Sample), y=reorder(NormalizedValue, Sample), group = Measure, color = Measure)) +
  #geom_point() +
  coord_polar(theta="x") +
  geom_polygon(fill = NA, linewidth = 0.8) + theme_minimal()+
  #theme(legend.position = "none")  +
  theme(axis.text.x = element_text(angle = 25, size = 7, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        plot.title = element_text(size=10, hjust = 0.5))+
  #ggtitle("SC: Quality score per each measure ")+
  xlim(SE_Quality$Sample)+
  scale_colour_manual(values=cbPalette_9)

SI_polygon <- ggplot(IT_Quality, aes(x = as.factor(Sample), y=reorder(NormalizedValue, Sample), group = Measure, color = Measure)) +
  #geom_point() +
  coord_polar(theta="x") +
  geom_polygon(fill = NA, linewidth = 0.8) + theme_minimal()+
  #theme(legend.position = "none")  +
  theme(axis.text.x = element_text(angle = 25, size = 7, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        plot.title = element_text(size=10, hjust = 0.5))+
  #ggtitle("SI: Quality score per each measure ")+
  xlim(IT_Quality$Sample)+
  scale_colour_manual(values=cbPalette)


#Arrange plots together
ggarrange(SC, SC_polygon,SI, SI_polygon,
          ncol = 2, nrow = 2)


#####PER SAMPLE##### Create the circular plot
IT_Quality$Adjusted_Score
Quality_overview <- normalized_data
normalized_data$NormalizedValue
IT <- ggplot(IT_Quality, aes(x = Sample, y=reorder(Adjusted_Score, Sample), group = Measure, color = Measure)) +
  #geom_point() +
  coord_polar(theta="x") +
  geom_polygon(fill = NA) + theme_minimal()+
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, size = 7, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        plot.title = element_text(size=10, hjust = 0.5))+
  ggtitle("IT20_007")+
  xlim(IT_Quality$Sample)+
  scale_colour_manual(values=cbPalette_9)
GR <- ggplot(GR_Quality, aes(x = Sample, y=reorder(Adjusted_Score, Sample), group = Measure, color = Measure)) +
  #geom_point() +
  coord_polar(theta="x") +
  geom_polygon(fill = NA) + theme_minimal()+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, size = 7, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        plot.title = element_text(size=10, hjust = 0.5))+
  ggtitle("GR01_002")+
  xlim(GR_Quality$Sample)+
  scale_colour_manual(values=cbPalette)
NO <- ggplot(NO_Quality, aes(x = Sample, y=reorder(Adjusted_Score, Sample), group = Measure, color = Measure)) +
  #geom_point() +
  coord_polar(theta="x") +
  geom_polygon(fill = NA) + theme_minimal()+
  #theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 30, size = 7, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        plot.title = element_text(size=10, hjust = 0.5))+
  ggtitle("NO01_002")+
  xlim(NO_Quality$Sample)+
  scale_colour_manual(values=cbPalette)
SE <- ggplot(SE_Quality, aes(x = as.factor(Sample), y=reorder(Adjusted_Score, Sample), group = Measure, color = Measure)) +
  #geom_point() +
  coord_polar(theta="x") +
  geom_polygon(fill = NA, linewidth = 0.8) + theme_minimal()+
  #theme(legend.position = "none")  +
  theme(axis.text.x = element_text(angle = 30, size = 7, hjust = 1),
        axis.title.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y =element_blank(),
        plot.title = element_text(size=10, hjust = 0.5))+
  ggtitle("SE02_002")+
  xlim(SE_Quality$Sample)+
  scale_colour_manual(values=cbPalette)


pdf("~/Studium/Köln_Biological Science/Master Thesis/Results/results.pdf")



