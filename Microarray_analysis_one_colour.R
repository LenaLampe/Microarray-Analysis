
library(statmod)
library(limma)
library(ggplot2)
library(utils)
library(reshape2)
library(annotate)
library(data.table)
library(plyr)
library(stats)
library(dplyr)


##############################################################################################
#Loading of miRNAs
##############################################################################################


setwd("file")

targets <- readTargets(file="targets.txt", row.names=NULL)

x <- read.maimages(targets, source="agilent",columns = list(G = "gProcessedSignal", 
                                                            Gb = "gProcessedBackground", Rb = "gIsPosAndSignif"),
                                                            green.only=TRUE)

y <- backgroundCorrect(x, method="normexp", offset=16)

y <- normalizeBetweenArrays(y, method="quantile")

y.ave <- avereps(y, ID=y$genes$ProbeName)


#####################################################################################################################
#Analysis for relative expression_Heatmap
#####################################################################################################################
# here I transform the Elist to one merged dataframe, which contains the log relative expression values of all samples
#####################################################################################################################


expression <- subset(y.ave$E)

gene_annotation <- subset(y.ave$genes)

expression <- as.data.frame(expression)

expression <- cbind(Row.Names = rownames(expression), expression)

colnames(expression)[1] <- "ProbeName"

row.names(expression) <- NULL

expression <- merge(expression, gene_annotation, by = "ProbeName")

write.table(expression, file="Relative_expression.txt", sep="\t", quote=FALSE)


############################################################################################################################
#Processing of the relative expression data table
############################################################################################################################


e.data <- read.delim("Relative_expression_all.txt")

e.data <- e.data[ , -which(names(e.data) %in% c("Row","Start","Col","ProbeUID","SystematicName","Description"))]


# I take the mean of all four replicates of the same sample type

e.data$Mean_Carcass_UI <- rowMeans(e.data[,2:5])

e.data$Mean_Head_UI <- rowMeans(e.data[,6:9])

e.data$Mean_Ovaries_UI <- rowMeans(e.data[,10:13])

e.data$Mean_Carcass_UF <- rowMeans(e.data[,14:17])

e.data$Mean_Head_UF <- rowMeans(e.data[,18:21])

e.data$Mean_Midgut_UF <- rowMeans(e.data[,22:25])

e.data$Mean_Ovaries_UF <- rowMeans(e.data[,26:29])

e.data$Mean_Carcass_I <- rowMeans(e.data[,30:33])

e.data$Mean_Head_I <- rowMeans(e.data[,34:37])

e.data$Mean_Ovaries_I <- rowMeans(e.data[,38:41])


#The data is saved in e.data table

write.table(e.data, file="Relative_expression_mean.txt", sep="\t", quote=FALSE)


# I delete all columns, which are not necessary for further analysis

e.data_Mean_only <- subset(e.data, select = c(1,42,44,45,46,47,48,49,50,51,52,53,54))


# I substract the negative Control from all columns (which are seperate microarrays) and I add 1 to not have to deal with negative values
e.data_Mean_only$Mean_Carcass_UI_minus_negative <- e.data_Mean_only$Mean_Carcass_UI - 4.772556 + 1
e.data_Mean_only$Mean_Head_UI_minus_negative <- e.data_Mean_only$Mean_Head_UI- 4.713432 +1
e.data_Mean_only$Mean_Ovaries_UI_minus_negative <- e.data_Mean_only$Mean_Ovaries_UI- 4.744095 +1
e.data_Mean_only$Mean_Carcass_UF_minus_negative <- e.data_Mean_only$Mean_Carcass_UF- 4.787483 +1
e.data_Mean_only$Mean_Head_UF_minus_negative <- e.data_Mean_only$Mean_Head_UF- 4.603362 +1
e.data_Mean_only$Mean_Midgut_UF_minus_negative <- e.data_Mean_only$Mean_Midgut_UF- 4.772428 +1
e.data_Mean_only$Mean_Ovaries_UF_minus_negative <- e.data_Mean_only$Mean_Ovaries_UF- 4.738908 +1 
e.data_Mean_only$Mean_Carcass_I_minus_negative <- e.data_Mean_only$Mean_Carcass_I- 4.738393 +1
e.data_Mean_only$Mean_Head_I_minus_negative <- e.data_Mean_only$Mean_Head_I- 4.750649 +1
e.data_Mean_only$Mean_Ovaries_I_minus_negative <- e.data_Mean_only$Mean_Ovaries_I- 4.829373 +1


# I save this data under a new file name

e.data_Mean_corrected_negative_control <- e.data_Mean_only

write.table(e.data_Mean_corrected_negative_control, file="Relative_expression_corrected_negative_control.txt", sep="\t", quote=FALSE)


# I delete the columns, which contain the values without the substraction of the negative control
e.data_Mean_corrected_negative_control <- subset(e.data_Mean_corrected_negative_control, select = c(1,2,3,14,15,16,17,18,19,20,21,22,23))


# I select all miRNA probes, which show in at least one sample a signal above the negative control
e.data_Mean_only_minus_negative_control_only_above_negative_control <- subset(e.data_Mean_corrected_negative_control, 
                                                                                Mean_Carcass_UI_minus_negative>1 |
                                                                                Mean_Head_UI_minus_negative>1 |
                                                                                Mean_Ovaries_UI_minus_negative>1 |
                                                                                Mean_Carcass_UF_minus_negative>1 |
                                                                                Mean_Head_UF_minus_negative> 1 |
                                                                                Mean_Midgut_UF_minus_negative>1 |
                                                                                Mean_Ovaries_UF_minus_negative>1 |
                                                                                Mean_Carcass_I_minus_negative>1 |
                                                                                Mean_Head_I_minus_negative>1 |
                                                                                Mean_Ovaries_I_minus_negative>1)

# I delete the row numbers
row.names(e.data_Mean_only_minus_negative_control_only_above_negative_control) <- NULL

write.table(e.data_Mean_only_minus_negative_control_only_above_negative_control, file="Relative_expression_all_above_negative_control.txt", sep="\t", quote=FALSE)


# I split the column of GeneName into the miRNA annotation and organism and then pick all the probes which are annotated as dipteran insect miRNAs
splitspecies <- colsplit(e.data_Mean_only_minus_negative_control_only_above_negative_control$GeneName, "-", c("Species", "miRNA"))

#I merge the two tables containing the organism information with the data
e.data_Mean_only_minus_negative_control_only_above_negative_control <- cbind(e.data_Mean_only_minus_negative_control_only_above_negative_control, splitspecies)

# I select all miRNA probes which relate to dipteran species
e.data_all_miRNAs <- subset(e.data_Mean_only_minus_negative_control_only_above_negative_control, Species == "aga"| Species == "aae"| Species =="ame"| Species == "cqu"| Species == "dme"| Species == "bmo"| Species == "tca")


write.table(e.data_all_miRNAs, file= "Relative_expression_insect_miRNAs_after_filtering_negative_control.txt", sep="\t", quote=FALSE)



###################### Now I blasted the sequences of all detected probes against all miRNAs detected and published on miRBase, Biryukova et al. and Castellano et. al.
######################This generates a list of probes that blasted against existing miRNAs including the number of nucleotides that the two sequences share. ##################################


positive_blasted_probes <- read.delim("Alignment_all_dipteran.csv", header = TRUE, sep = ',')

positive_blasted_probes <- subset(positive_blasted_probes, select = c(ProbeName, nucleotide_overlap, miRNA_annotation))

# Now I merge the e.data_all_miRNA table withe the positive_blasted_probes and I have a dataframe which contains the information how
# many probe nucleotides match the miRNA

miRNA_blasted <- merge(e.data_all_miRNAs, positive_blasted_probes, by = "ProbeName", all.x = TRUE)

# I now pick all probes, which have at least 14 nucleotides overlap with the miRNA

miRNA_blasted_cut_off_14 <- subset(miRNA_blasted, nucleotide_overlap >= 14)

write.table(miRNA_blasted_cut_off_14, file="Relative_expression_filtered.txt", sep = '\t', quote = FALSE)


sapply(miRNA_blasted_cut_off_14, class)


final_miRNA_set <- ddply(miRNA_blasted_cut_off_14,"miRNA_annotation",numcolwise(mean))

###Now I subset only the samples which were unfed


final_miRNA_set_Unfed <- subset(final_miRNA_set, select = c (miRNA_annotation,  
                                                             Mean_Carcass_UF_minus_negative, Mean_Head_UF_minus_negative, 
                                                             Mean_Midgut_UF_minus_negative, Mean_Ovaries_UF_minus_negative))

final_miRNA_set_Unfed <- plyr::rename(final_miRNA_set_Unfed, c("miRNA_annotation" = "miRNA", "Mean_Carcass_UF_minus_negative" = "Carcass", "Mean_Head_UF_minus_negative" = "Head",
                                                               "Mean_Midgut_UF_minus_negative" = "Midgut", "Mean_Ovaries_UF_minus_negative" = "Ovary"))

write.table(final_miRNA_set_Unfed, file="final_miRNA_set_Unfed.txt", sep="\t", quote=FALSE)

#####################################################################################################################################
#Heatmap for Unfed 



final_miRNA_set_Unfed <- read.delim("final_miRNA_set_Unfed.csv", header = TRUE, sep = ';', dec = ',')

splitspecies2 <- colsplit(final_miRNA_set_Unfed$miRNA, "-", c("Species", "MiRNA"))

final_miRNA_set_Unfed <- cbind(final_miRNA_set_Unfed, splitspecies2)

final_miRNA_set_Unfed$miRNA <- NULL
final_miRNA_set_Unfed$Species <- NULL

final_miRNA_set_Unfed <- final_miRNA_set_Unfed[,c(5,2,3,1,4)]

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

library("plyr")
library("reshape2")
library("heatmap.plus")

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


data <- final_miRNA_set_Unfed
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names 

#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette

my_palette <- colorRampPalette(c("snow","dodgerblue", "sandybrown", "saddlebrown"))(n = 499)


#(optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,1.4,length=100),               # for red
               seq(1.5,2.8,length=100),             # for yellow
               seq(2.9,4.2,length=100),
               seq(4.3,5.6,length=100),
               seq(5.7,7,length=100))               # for green


pdf ("UF_heatmap.pdf",   # c
     width = 6,        # 5 x 300 pixels
     height = 15)        # smaller font size


heatmap.2(mat_data,        
          key= TRUE,
          keysize = 2,
          key.xlab = "relative expression",     
          density.info="none",                
          trace="none",                       
          margins =c(12,12), 
          cexCol = 1.0 , 
          cexRow = 0.45,
          col=my_palette,                      
          breaks=col_breaks,                  
          hclustfun = hclust, 
          dendrogram="both",
          offsetRow = 0,
          srtCol = 45,
          Colv=T)                


dev.off()                       # close the PNG device


#########################################################################################################################
#######Heatmap for Blood fed mosquitoes



#########################################################
### B) Reading in data and transform it into matrix format
#########################################################


###Now I subset only the samples which were blood fed


final_miRNA_set_BM <- subset(final_miRNA_set, select = c (miRNA_annotation,  
                                                          Mean_Carcass_UI_minus_negative, Mean_Head_UI_minus_negative, 
                                                          Mean_Ovaries_UI_minus_negative))

final_miRNA_set_BM <- plyr::rename(final_miRNA_set_BM, c("miRNA_annotation" = "miRNA","Mean_Carcass_UI_minus_negative" = "Carcass_BM", "Mean_Head_UI_minus_negative" = "Head_BM",
                                                         "Mean_Ovaries_UI_minus_negative" = "Ovary_BM"))

splitspecies3 <- colsplit(final_miRNA_set_BM$miRNA, "-", c("Species", "MiRNA"))

final_miRNA_set_BM <- cbind(final_miRNA_set_BM, splitspecies3)

final_miRNA_set_BM$miRNA <- NULL
final_miRNA_set_BM$Species <- NULL

final_miRNA_set_BM <- final_miRNA_set_BM[,c(4,2,1,3)]

write.table(final_miRNA_set_BM, file="final_miRNA_set_BM.txt", sep="\t", quote=FALSE)

data <- final_miRNA_set_BM
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names 

#########################################################
### C) Customizing and plotting the heat map
#########################################################


#creates a own color palette from red to green
my_palette <- colorRampPalette(c("snow","lightskyblue", "royalblue", "darkblue"))(n = 499)


#(optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,1.4,length=100),               # for red
               seq(1.5,2.8,length=100),             # for yellow
               seq(2.9,4.2,length=100),
               seq(4.3,5.6,length=100),
               seq(5.7,7,length=100))               # for green

#creates a 5 x 5 inch image
pdf ("BM_heatmap.pdf",   # c
     width = 6,        # 5 x 300 pixels
     height = 12)        # smaller font size

heatmap.2(mat_data, 
          #main = "Unfed",        # heat map title
          key= TRUE,             # change font color of cell labels to black
          density.info="none",   # turns off density plot inside color legend
          trace="none",          # turns off trace lines inside the heat map
          margins =c(12,12), 
          cexCol = 1.5 , 
          cexRow = 0.45,
          col=my_palette,        # use on color palette defined earlier 
          breaks=col_breaks,     # enable color transition at specified limits 
          hclustfun = hclust, 
          dendrogram="row",
          offsetRow = 0,
          srtCol = 45,
          Colv=F)                # turn off column clustering


dev.off()                        # close the PNG device



#####################################################################################################################
#Differential miRNA expression by blood feeding and Infection
#####################################################################################################################


Treat <- factor(paste(targets$Condition,targets$Tissue,sep="."))

design <- model.matrix(~0+Treat)

colnames(design) <- levels(Treat)

corfit <- duplicateCorrelation(y.ave,design,block=targets$Subject)

corfit$consensus

fit <- lmFit(y.ave,design,block=targets$Subject, correlation = corfit$consenus)


contrast.matrix <- makeContrasts("Uninfected.Carcass-Unfed.Carcass", "Infected.Carcass-Uninfected.Carcass",
                                 "Uninfected.Ovaries-Unfed.Ovaries", "Infected.Ovary-Uninfected.Ovaries",
                                 "Uninfected.Head-Unfed.Head", "Infected.Head-Uninfected.Head",
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


#########################################################################
#Output for differential expression
#########################################################################



Blood_meal_Carcass <- topTable(fit2, adjust = "BH",coef="Uninfected.Carcass-Unfed.Carcass", genelist=y.ave$genes, number=Inf)
write.table(Blood_meal_Carcass, file="C:/R data/data/Analysis_final_09_09_16/Differential_expression_BM/Carcass/BM_Carcass.txt", sep="\t", quote=FALSE)

Infection_Carcass <- topTable(fit2, adjust="BH", coef="Infected.Carcass-Uninfected.Carcass", genelist=y.ave$genes, number=Inf)
write.table(Infection_Carcass, file="C:/R data/data/Analysis_final_09_09_16/Differential_expression_I/Carcass/I_Carcass.txt", sep="\t", quote=FALSE)



Blood_meal_Head <- topTable(fit2, adjust="BH", coef="Uninfected.Head-Unfed.Head", genelist=y.ave$genes, number=Inf)
write.table(Blood_meal_Head, file="C:/R data/data/Analysis_final_09_09_16/Differential_expression_BM/Head/BM_Head.txt", sep="\t", quote=FALSE)

Infection_Head <- topTable(fit2, adjust="BH", coef="Infected.Head-Uninfected.Head", genelist=y.ave$genes, number=Inf)
write.table(Infection_Head, file="C:/R data/data/Analysis_final_09_09_16/Differential_expression_I/Head/I_Head.txt", sep="\t", quote=FALSE)



Blood_meal_Ovary <- topTable(fit2, adjust="holm", coef="Uninfected.Ovaries-Unfed.Ovaries", genelist=y.ave$genes, number=Inf)
write.table(Blood_meal_Ovary, file="C:/R data/data/Analysis_final_09_09_16/Differential_expression_BM/Ovary/BM_Ovary.txt", sep="\t", quote=FALSE)

Infection_Ovary <- topTable(fit2, adjust="BH", coef="Infected.Ovary-Uninfected.Ovaries", genelist=y.ave$genes, number=Inf)
write.table(Infection_Ovary, file="C:/R data/data/Analysis_final_09_09_16/Differential_expression_I/Ovary/I_Ovary.txt", sep="\t", quote=FALSE)


##############################################################################################################################################
#########Volcano Plots for Blood meals
##########################################################################################################################################



#########Differential expression of miRNAs in the Carcass after blood meal intake, including only the miRNA probes which match at least 14 nucleotides the miRNA

setwd("C:/R data/data/Analysis_final_09_09_16/Differential_expression_BM/Carcass")

Probe_ID_included <- subset(miRNA_blasted_cut_off_14, select = c(ProbeName, miRNA_annotation, nucleotide_overlap))

Blood_meal_Carcass <- merge(Blood_meal_Carcass, Probe_ID_included, by = "ProbeName")


##Highlight genes 
Blood_meal_Carcass$threshold = as.factor(abs(Blood_meal_Carcass$logFC) > 0.5849 & Blood_meal_Carcass$adj.P.Val < 0.05)

##Construct the plot object
Blood_meal_Carcass_vplot = ggplot(data=Blood_meal_Carcass, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(data = subset(Blood_meal_Carcass, threshold == TRUE), alpha=0.8, size=5, colour = "red") +
  geom_point(data = subset(Blood_meal_Carcass, threshold == FALSE), alpha=0.8, size=4, colour = "grey50")+
  geom_text(data = subset(Blood_meal_Carcass, threshold == TRUE), aes(label = miRNA_annotation ), size = 4, hjust = 1.2, vjust = 1.1) +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.ticks = element_line (colour = "black",size = 1,),  
        axis.title.x = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.title.y = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.text  = element_text( colour = "black", vjust=0.5, size=18)) +
  scale_x_continuous(limit = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(limit = c(0, 8), breaks = c(0,1,2,3,4,5,6,7,8))+
  xlab("log2 (Fold Change)") + ylab("-log10 (p-value)") +
  geom_hline(yintercept=1.3, colour = "red") + geom_vline(xintercept=0.5849, colour = "red") + geom_vline(xintercept=-0.5849, colour = "red")
Blood_meal_Carcass_vplot
dev.off

ggsave(filename = "Blood_meal_Carcass.pdf", plot = Blood_meal_Carcass_vplot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c( "cm"),
       dpi = 300, limitsize = TRUE)


############### Volcano plot for Head_BM

setwd("C:/R data/data/Analysis_final_09_09_16/Differential_expression_BM/Head")


Blood_meal_Head <- merge(Blood_meal_Head, Probe_ID_included, by = "ProbeName")

##Highlight genes 

Blood_meal_Head$threshold = as.factor(abs(Blood_meal_Head$logFC) > 0.5849 & Blood_meal_Head$adj.P.Val < 0.05)


##Construct the plot object
Blood_meal_Head_vplot = ggplot(data=Blood_meal_Head, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(data = subset(Blood_meal_Head, threshold == TRUE), alpha=0.8, size=5, colour = "red") +
  geom_point(data = subset(Blood_meal_Head, threshold == FALSE), alpha=0.8, size=4, colour = "grey50")+
  geom_text(data = subset(Blood_meal_Head, threshold == TRUE), aes(label = miRNA_annotation ), size = 4, hjust = 1.2, vjust = 1.1) +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.ticks = element_line (size = 1, colour = "black"),  
        axis.title.x = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.title.y = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.text  = element_text( colour = "black", vjust=0.5, size=18)) +
  scale_x_continuous(limit = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(limit = c(0, 8), breaks = c(0,1,2,3,4,5,6,7,8))+
  xlab("log2 (Fold Change)") + ylab("-log10 (p-value)") +
  geom_hline(yintercept=1.3, colour = "red") + geom_vline(xintercept=0.5849, colour = "red") + geom_vline(xintercept=-0.5849, colour = "red")
Blood_meal_Head_vplot

ggsave(filename = "Blood_meal_Head.pdf", plot = Blood_meal_Head_vplot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c( "cm"),
       dpi = 300, limitsize = TRUE)


#########Volcano plot for Ovary_BM

##Highlight genes 

setwd("C:/R data/data/Analysis_final_09_09_16/Differential_expression_BM/Ovary")

Blood_meal_Ovary <- merge(Blood_meal_Ovary, Probe_ID_included, by = "ProbeName")

Blood_meal_Ovary$threshold = as.factor(abs(Blood_meal_Ovary$logFC) > 0.5849 & Blood_meal_Ovary$adj.P.Val < 0.05)

##Construct the plot object
Blood_meal_Ovary_vplot = ggplot(data=Blood_meal_Ovary, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(data = subset(Blood_meal_Ovary, threshold == TRUE), alpha=0.8, size=5, colour = "red") +
  geom_point(data = subset(Blood_meal_Ovary, threshold == FALSE), alpha=0.8, size=4, colour = "grey50")+
  geom_text(data = subset(Blood_meal_Ovary, threshold == TRUE), aes(label = miRNA_annotation ), size = 4, hjust = 1.2, vjust = 1.1) +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.ticks = element_line (size = 1, colour = "black"),  
        axis.title.x = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.title.y = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.text  = element_text( colour = "black", vjust=0.5, size=18)) +
  scale_x_continuous(limit = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(limit = c(0, 8), breaks = c(0,1,2,3,4,5,6,7,8))+
  xlab("log2 (Fold Change)") + ylab("-log10 (p-value)") +
  geom_hline(yintercept=1.3, colour = "red") + geom_vline(xintercept=0.5849, colour = "red") + geom_vline(xintercept=-0.5849, colour = "red")
Blood_meal_Ovary_vplot

ggsave(filename = "Blood_meal_Ovary.pdf", plot = Blood_meal_Ovary_vplot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c( "cm"),
       dpi = 300, limitsize = TRUE)



##################################################################################
#################Volcano Plot for Infection
###################################################################################

setwd("C:/R data/data/Analysis_final_09_09_16/Differential_expression_I/Carcass")


Probe_ID_included <- subset(miRNA_blasted_cut_off_14, select = c(ProbeName, miRNA_annotation, nucleotide_overlap))


Infection_Carcass <- merge(Infection_Carcass, Probe_ID_included, by = "ProbeName")


##Highlight genes 
Infection_Carcass$threshold = as.factor(abs(Infection_Carcass$logFC) > 0.5849 & Infection_Carcass$adj.P.Val < 0.05)

##Construct the plot object
Infection_Carcass_vplot = ggplot(data=Infection_Carcass, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(data = subset(Infection_Carcass, threshold == TRUE), alpha=0.8, size=5, colour = "red") +
  geom_point(data = subset(Infection_Carcass, threshold == FALSE), alpha=0.8, size=4, colour = "grey50")+
  geom_text(data = subset(Infection_Carcass, threshold == TRUE), aes(label = miRNA_annotation ), size = 4, hjust = 1.2, vjust = 1.1) +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.ticks = element_line (colour = "black",size = 1,),  
        axis.title.x = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.title.y = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.text  = element_text( colour = "black", vjust=0.5, size=18)) +
  scale_x_continuous(limit = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(limit = c(0, 8), breaks = c(0,1,2,3,4,5,6,7,8))+
  xlab("log2 (Fold Change)") + ylab("-log10 (p-value)") +
  geom_hline(yintercept=1.3, colour = "red") + geom_vline(xintercept=0.5849, colour = "red") + geom_vline(xintercept=-0.5849, colour = "red")
Infection_Carcass_vplot
dev.off

ggsave(filename = "Infection_Carcass.pdf", plot = Infection_Carcass_vplot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c( "cm"),
       dpi = 300, limitsize = TRUE)

############### Volcano plot for Head_BM

setwd("C:/R data/data/Analysis_final_09_09_16/Differential_expression_I/Head")

Infection_Head <- merge(Infection_Head, Probe_ID_included, by = "ProbeName")

##Highlight genes 

Infection_Head$threshold = as.factor(abs(Infection_Head$logFC) > 0.5849 & Infection_Head$adj.P.Val < 0.05)


##Construct the plot object
Infection_Head_vplot = ggplot(data=Infection_Head, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(data = subset(Infection_Head, threshold == TRUE), alpha=0.8, size=5, colour = "red") +
  geom_point(data = subset(Infection_Head, threshold == FALSE), alpha=0.8, size=4, colour = "grey50")+
  geom_text(data = subset(Infection_Head, threshold == TRUE), aes(label = miRNA_annotation ), size = 4, hjust = 1.2, vjust = 1.1) +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.ticks = element_line (colour = "black",size = 1,),  
        axis.title.x = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.title.y = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.text  = element_text( colour = "black", vjust=0.5, size=18)) +
  scale_x_continuous(limit = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(limit = c(0, 8), breaks = c(0,1,2,3,4,5,6,7,8))+
  xlab("log2 (Fold Change)") + ylab("-log10 (p-value)") +
  geom_hline(yintercept=1.3, colour = "red") + geom_vline(xintercept=0.5849, colour = "red") + geom_vline(xintercept=-0.5849, colour = "red")
Infection_Head_vplot
dev.off

ggsave(filename = "Infection_Head.pdf", plot = Infection_Head_vplot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c( "cm"),
       dpi = 300, limitsize = TRUE)
#########Volcano plot for Ovary_BM

##Highlight genes 

setwd("C:/R data/data/Analysis_final_09_09_16/Differential_expression_I/Ovary")

Infection_Ovary <- merge(Infection_Ovary, Probe_ID_included, by = "ProbeName")

Infection_Ovary$threshold = as.factor(abs(Infection_Ovary$logFC) > 0.5849 & Infection_Ovary$adj.P.Val < 0.05)

##Construct the plot object
Infection_Ovary_vplot = ggplot(data=Infection_Ovary, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(data = subset(Infection_Ovary, threshold == TRUE), alpha=0.8, size=5, colour = "red") +
  geom_point(data = subset(Infection_Ovary, threshold == FALSE), alpha=0.8, size=4, colour = "grey50")+
  geom_text(data = subset(Infection_Ovary, threshold == TRUE), aes(label = miRNA_annotation ), size = 4, hjust = 1.2, vjust = 1.1) +
  theme(legend.position = "right", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.ticks = element_line (colour = "black",size = 1,),  
        axis.title.x = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.title.y = element_text( size= 18, hjust = 0.5, vjust = 1.5),
        axis.text  = element_text( colour = "black", vjust=0.5, size=18)) +
  scale_x_continuous(limit = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(limit = c(0, 8), breaks = c(0,1,2,3,4,5,6,7,8))+
  xlab("log2 (Fold Change)") + ylab("-log10 (p-value)") +
  geom_hline(yintercept=1.3, colour = "red") + geom_vline(xintercept=0.5849, colour = "red") + geom_vline(xintercept=-0.5849, colour = "red")
Infection_Ovary_vplot
dev.off

ggsave(filename = "Infection_Ovary.pdf", plot = Infection_Ovary_vplot, device = NULL, path = NULL,
       scale = 1, width = 12, height = 12, units = c( "cm"),
       dpi = 300, limitsize = TRUE)
