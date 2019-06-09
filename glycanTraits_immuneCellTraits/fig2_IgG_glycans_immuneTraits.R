## Tiphaine Martin - 6/6/2019
## Script developped for the paper doi: https://doi.org/10.1101/662957

## if you use this script as template, please cite this paper
## doi: https://doi.org/10.1101/662957

############################# 
#        Figure 2
#############################
##Correlation between immune traits associated with IgG
library(ComplexHeatmap)
require(circlize)
igG_immune <-read.table(file="IgG_ImmuneTraits2.csv",header=TRUE,sep=",")
rownames(igG_immune) <- igG_immune$X
igG_immune<- igG_immune[,-1]
igG_immune <- as.matrix(igG_immune)
coul=col3(100)
#coul = colorRampPalette(brewer.pal(8, "RdBu"))(25)
heatmap(igG_immune,col= coul)
colRwB=colorRamp2(c(-0.5, -0.4,-0.3,-0.2,-0.1,
                    0, 
                    0.1,0.2,0.3,0.4,0.5), 
                  c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                    "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                    "#D6604D", "#B2182B", "#67001F"))

Heatmap(igG_immune,col=colRwB,color_space = colRwB)

data        = data.matrix(igG_immune)
distance    = dist(data)
cluster     = hclust(distance, method="ward.D2")
dendrogram  = as.dendrogram(cluster)
Rowv        = rowMeans(data, na.rm = T)
dendrogram  = reorder(dendrogram, Rowv)
rowInd = order.dendrogram(dendrogram)
rowInd_names <- rownames(igG_immune)[rowInd]
Heatmap(igG_immune,col=colRwB,color_space = colRwB, cluster_rows=dendrogram)

#### Correlation 17 IgG N-glycan structures
library("RColorBrewer")
library("corrplot")

mycorrelation <- read.csv(file="/Users/tiphainecmartin/Documents/Projects_KCL/AITD_glycosylation/Ewa/coMET/ST3_AITD_antiTPO_GlycomicWAS_1/glycans_poland_igg_global_combat_scale_cleaned2018.csv",header=TRUE)

col3 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))  

dim(mycorrelation)
subIgG <- colnames(mycorrelation)[which(colnames(mycorrelation) %in% c('IGP42','IGP21','IGP2','IGP36','IGP45','IGP46',
                                                                       'IGP59','IGP48','IGP8','IGP56','IGP58','IGP15',
                                                                       'IGP60','IGP7','IGP33','IGP62','IGP63'))]

submycorrelation <- mycorrelation[,subIgG]
dim(submycorrelation)

submycorrelation <- submycorrelation[, c('IGP2', 'IGP7', 'IGP8','IGP15','IGP21','IGP33',
                                         'IGP36','IGP42','IGP45','IGP46','IGP48',
                                         'IGP56','IGP58','IGP59',
                                         'IGP60','IGP62','IGP63')]

Msub <-cor(submycorrelation)
res1sub <- cor.mtest(submycorrelation, conf.level = .95)
corrplot(Msub, type = "upper", method = "ellipse",tl.cex =0.75, tl.col = "black",
         p.mat = res1sub$p, insig = "blank", sig.level = .05, col=col3(200),order = "hclust")


### correlation between 51 immune cell traits
data_path="/figures_papers/data/"
#data clean
load(paste0(data_path,"ds.age.and.batch.corrected.quantile.normalised_all.RData"))
dim(ds.age.and.batch.corrected.quantile.normalised_interesting_date_clean)

immune_features <- t(ds.age.and.batch.corrected.quantile.normalised_interesting_date_clean[,-c(1:6)])
dim(immune_features)

immunes51 <-read.table(file="/figures_papers/data/list_51immuneTraits.csv",header=TRUE,sep=",")
dim(immunes51)
immunes51_order <- immunes51[match(rowInd_names, immunes51$Subset.name),]
dim(immunes51_order)

immune_features_sub51 <- immune_features[which(rownames(immune_features)%in% immunes51_order$immuneTrait.ID),]
dim(immune_features_sub51)
immune_features_sub51_order <- immune_features_sub51[match(immunes51_order$immuneTrait.ID, rownames(immune_features_sub51)),]
rownames(immune_features_sub51_order)
rownames(immune_features_sub51_order) <- immunes51_order$Subset.name
rownames(immune_features_sub51_order)
rownames(immune_features_sub51_order) <- gsub("^16\\+56\\/","",rownames(immune_features_sub51_order))

immune_features_sub51_v1 <- immune_features_sub51_order[,-which(colnames(immune_features_sub51_order) %in% c(159,197,260,468,499))]
dim(immune_features_sub51_v1)

timmune_features_sub51_v1 <- t(immune_features_sub51_v1)
dim(timmune_features_sub51_v1)
immune_features_sub51_NA <- timmune_features_sub51_v1[complete.cases(timmune_features_sub51_v1), ]
dim(immune_features_sub51_NA)

MsubImm <-cor(immune_features_sub51_NA)
dim(MsubImm)
res1subImm <- cor.mtest(immune_features_sub51_NA, conf.level = .95)
corrplot(MsubImm, type = "upper", method = "ellipse",tl.cex =0.75, tl.col = "black",
         p.mat = res1subImm$p, insig = "blank", sig.level = .05, col=col3(200))
