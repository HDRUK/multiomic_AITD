## Tiphaine Martin - 6/6/2019
## Script developped for the paper doi: https://doi.org/10.1101/662957

## if you use this script as template, please cite this paper
## doi: https://doi.org/10.1101/662957

############################# 
#        Figure 3
#############################
## correlation between immune cell traits - UMAP plot
use_python("/anaconda3/bin/python3.6")
devtools::install_github("ropenscilabs/umapr")
library(umapr)
library(tidyverse)

data_path="/figures_papers/data/"
#data clean
load(paste0(data_path,"ds.age.and.batch.corrected.quantile.normalised_all.RData"))
dim(ds.age.and.batch.corrected.quantile.normalised_interesting_date_clean)

immune_features <- t(ds.age.and.batch.corrected.quantile.normalised_interesting_date_clean[,-c(1:6)])
dim(immune_features)
immune_features_v1 <- immune_features[,-which(colnames(immune_features) %in% c(159,197,260,468,499))]
dim(immune_features_v1)

## replace NA by 0
immune_features_v1[is.na(immune_features_v1)] <- 0
dim(immune_features_v1)

# run UMAP algorithm
pembedding_immune <- umap(immune_features_v1)

#add annotation
immune <-read.table(file="/figures_papers/data/immune_traits_annotations.csv",header=TRUE,sep=",")
pembedding_immune$immuneTrait <- rownames(pembedding_immune)
pembedding_immune_annot <- merge(pembedding_immune,immune, by="immuneTrait")
pembedding_immune_annot$annotation <- pembedding_immune_annot$Lineage
pembedding_immune_annot$annotation <- as.character(pembedding_immune_annot$annotation)
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "NA")),"annotation"] <- "other"
pembedding_immune_annot[which(is.na(pembedding_immune_annot$annotation)),"annotation"] <- "other"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "%")),"annotation"] <- "other"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "08")),"annotation"] <- "other"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "16+56-")),"annotation"] <- "NK:16+"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "16+56")),"annotation"] <- "NK:16+56+"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "11c+")),"annotation"] <- "11c+"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "14-")),"annotation"] <- "14-"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "Vg9+")),"annotation"] <- "Vg9+/Vd2"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation, "4+")),"annotation"] <- "4+"

pembedding_immune_annot$annotation2  <- pembedding_immune_annot$annotation 
pembedding_immune_annot$annotation2 <- as.character(pembedding_immune_annot$annotation2)
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "B cells")),"annotation2"] <- "B cells"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "CD4")),"annotation2"] <- "CD4"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "4")),"annotation2"] <- "CD4"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "56")),"annotation2"] <- "NK"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "NK")),"annotation2"] <- "NK"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "1")),"annotation2"] <- "other"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "CD3")),"annotation2"] <- "CD3"
pembedding_immune_annot[which(startsWith(pembedding_immune_annot$annotation2, "Vg")),"annotation2"] <- "other"


# plot the result
ggplot(pembedding_immune, aes(UMAP1, UMAP2)) + geom_point() +theme_bw()

pembedding_immune_annot$annotation <- as.factor(pembedding_immune_annot$annotation)
ggplot(pembedding_immune_annot,aes(UMAP1, UMAP2, color = annotation)) + 
  geom_point() +theme_bw()

pembedding_immune_annot$annotation2 <- as.factor(pembedding_immune_annot$annotation2)
ggplot(pembedding_immune_annot,aes(UMAP1, UMAP2, color = annotation2)) + 
  geom_point() +theme_bw()

## coMET plot
library(coMET)
library("rtracklayer")

chrom <- "chr6" 
start <- 31215228
end <- 31491679 # 31541779
gen <- "hg19" 
strand <- "*"
BROWSER.SESSION="UCSC" 
mySession <- browserSession(BROWSER.SESSION) 
genome(mySession) <- gen

itrack <- IdeogramTrack(genome = gen, chromosome = chrom)
gtrack <- GenomeAxisTrack()

ENSEMBLtrack <-genes_ENSEMBL(gen,chrom,start,end,showId=TRUE)
#HLA-C=ENSG00000204525
geneTrackShow_HLAC <- ENSEMBLtrack[gene(ENSEMBLtrack) %in% c("ENSG00000204525")]
start_exonHLAC <- start(geneTrackShow_HLAC@range@ranges)
end_exonHLAC <- end(geneTrackShow_HLAC@range@ranges)

#MIC-A=ENSG00000204520
geneTrackShow_MICA <- ENSEMBLtrack[gene(ENSEMBLtrack) %in% c("ENSG00000204520")]
start_exonMICA <- start(geneTrackShow_MICA@range@ranges)
end_exonMICA <- end(geneTrackShow_MICA@range@ranges)

#MIC-B=ENSG00000204516
geneTrackShow_MICB <- ENSEMBLtrack[gene(ENSEMBLtrack) %in% c("ENSG00000204516")]
start_exonMICB <- start(geneTrackShow_MICB@range@ranges)
end_exonMICB <- end(geneTrackShow_MICB@range@ranges)

#HCP5=ENSG00000206337
geneTrackShow_HCP5 <- ENSEMBLtrack[gene(ENSEMBLtrack) %in% c("ENSG00000206337")]
start_exonHCP5 <- start(geneTrackShow_HCP5@range@ranges)
end_exonHCP5 <- end(geneTrackShow_HCP5@range@ranges)

interestfeatures <- rbind(cbind(start_exonHLAC,end_exonHLAC,"HLA-C"),
                          cbind(start_exonMICA,end_exonMICA,"MIC-A"),
                          cbind(start_exonMICB,end_exonMICB,"MIC-B"),
                          cbind(start_exonHCP5,end_exonHCP5,"HCP5"))
interestcolor <- list("HLA-C"="red","MIC-A"="red","MIC-B"="red",
                      "HCP5"="red")

genetrack_highlighted <-interestGenes_ENSEMBL(gen,chrom,start,end,
                                              interestfeatures, interestcolor,
                                              showId=TRUE)


track.name="Broad ChromHMM" 
tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name)) 
table.name<-tablestrack[1] 
nameTissue<-"GM12878"
chromhmmPattern<-chromatinHMMAll_UCSC(gen,chrom,start,end,mySession,color='coMET',
                                      nameTissue) 

##DNASE
dnasetrack<-DNAse_UCSC(gen,chrom,start,end,mySession)
##CpG island
cpgIstrack <- cpgIslands_UCSC(gen,chrom,start,end)

##SNP
getOption("Gviz.scheme")
scheme <- getScheme()
scheme$AnnotationTrack$size <- 2

dirFiles <- "/Users/tiphainecmartin/Documents/Projects_KCL/AITD_Immunefeatures/version5/tables_figures_suppl_ter/plot/data_fig1a/"
AITD_snp_file <- paste0(dirFiles,"AITD_SNPs.txt")
read.table(AITD_snp_file,header = TRUE, sep = "\t")
AITDsnp <- eQTL(gen,chrom,start, end, AITD_snp_file, 
                featureDisplay = "SNP" )

displayPars(AITDsnp) <- list(size = 2) 

immune_snp_file <- paste0(dirFiles,"ImmuneTraits_SNPs.txt")
read.table(immune_snp_file,header = TRUE, sep = "\t")
Immunesnp <- eQTL(gen,chrom,start, end, immune_snp_file, 
                  featureDisplay = "SNP" )
displayPars(Immunesnp) <- list(size = 2)

IgG_snp_file <- paste0(dirFiles,"IgG_SNPs.txt")
read.table(IgG_snp_file,header = TRUE, sep = "\t")
IgGsnp <- eQTL(gen,chrom,start, end, IgG_snp_file, 
                  featureDisplay = "SNP" )
displayPars(IgGsnp) <- list(size = 2)

## eQTL
eQTL_thyroid_file <- paste0(dirFiles,"eQTL_thyroid.txt")
read.table(eQTL_thyroid_file,header = TRUE, sep = "\t")
eQTL_thyroid <- eQTL(gen,chrom,start, end, eQTL_thyroid_file, 
                  featureDisplay = "all", 
                  type_stacking="squish",just_group="above")
displayPars(eQTL_thyroid) <- list(size = 2)

eQTL_blood_file <- paste0(dirFiles,"eQTL_blood.txt")
read.table(eQTL_blood_file,header = TRUE, sep = "\t")
eQTL_blood <- eQTL(gen,chrom,start, end, eQTL_blood_file, 
                     featureDisplay = "all", 
                     type_stacking="squish",just_group="above")
displayPars(eQTL_blood) <- list(size = 2)

### chromHMM Thyroid
options(Gviz.scheme = "default")
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_40_19_segments.bed > CEMT_40_19_segments_AITD.bed
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_41_19_segments.bed > CEMT_41_19_segments_AITD.bed
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_42_19_segments.bed > CEMT_42_19_segments_AITD.bed
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_43_19_segments.bed > CEMT_43_19_segments_AITD.bed
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_44_19_segments.bed > CEMT_44_19_segments_AITD.bed
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_45_19_segments.bed > CEMT_45_19_segments_AITD.bed
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_86_19_segments.bed > CEMT_86_19_segments_AITD.bed
#awk -F"\t" '$1 == "chr6" && $2 > 31215228 && $3 < 31501779 {print $0}' CEMT_87_19_segments.bed > CEMT_87_19_segments_AITD.bed

colorchromHMM <- list('U1' = '#FF0000', 'U2' = '#FF4500','U3' = '#FF4500', 
                      'U4' = '#FF4500','U5' = '#008000','U6' = '#008000',
                      'U7' = '#006400', 'U8' = '#C2FF05','U9' = '#C2FF05',
                      'U10' = '#FFC34D','U11' = '#FFFF00', 'U12' = '#FFFF00',
                      'U13' = '#66CDAA','U14' = '#8A91D0','U15' = '#8A91D0',
                      'U16' = '#BDB76B', 'U17' = '#C0C0C0', 'U18' = '#808080',
                      'U19' = '#F7F7F7')

# data 40
chromHMMTrack_40_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_40_19_segments_AITD.bed")
chromHMMTrack_40_tab <- read.table(file=chromHMMTrack_40_path,header=FALSE,sep="\t")

chromHMMTrack_40_mat <- data.frame(chr="chr6",start=chromHMMTrack_40_tab[,2],
                                   end=chromHMMTrack_40_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_40_tab[,4],
                                   id=1:nrow(chromHMMTrack_40_tab),
                                   group=1:nrow(chromHMMTrack_40_tab))
chromHMMTrack_40_gr <- makeGRangesFromDataFrame(chromHMMTrack_40_mat, TRUE)
chromHMMTrack_40 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_40_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT40 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_40) <- colorchromHMM
displayPars(chromHMMTrack_40) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 
# data 41
chromHMMTrack_41_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_41_19_segments_AITD.bed")
chromHMMTrack_41_tab <- read.table(file=chromHMMTrack_41_path,header=FALSE,sep="\t")

chromHMMTrack_41_mat <- data.frame(chr="chr6",start=chromHMMTrack_41_tab[,2],
                                   end=chromHMMTrack_41_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_41_tab[,4],
                                   id=1:nrow(chromHMMTrack_41_tab),
                                   group=1:nrow(chromHMMTrack_41_tab))
chromHMMTrack_41_gr <- makeGRangesFromDataFrame(chromHMMTrack_41_mat, TRUE)
chromHMMTrack_41 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_41_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT41 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_41) <- colorchromHMM
displayPars(chromHMMTrack_41) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 

# data 42
chromHMMTrack_42_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_42_19_segments_AITD.bed")
chromHMMTrack_42_tab <- read.table(file=chromHMMTrack_42_path,header=FALSE,sep="\t")

chromHMMTrack_42_mat <- data.frame(chr="chr6",start=chromHMMTrack_42_tab[,2],
                                   end=chromHMMTrack_42_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_42_tab[,4],
                                   id=1:nrow(chromHMMTrack_42_tab),
                                   group=1:nrow(chromHMMTrack_42_tab))
chromHMMTrack_42_gr <- makeGRangesFromDataFrame(chromHMMTrack_42_mat, TRUE)
chromHMMTrack_42 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_42_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT42 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_42) <- colorchromHMM
displayPars(chromHMMTrack_42) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 

# data 43
chromHMMTrack_43_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_43_19_segments_AITD.bed")
chromHMMTrack_43_tab <- read.table(file=chromHMMTrack_43_path,header=FALSE,sep="\t")

chromHMMTrack_43_mat <- data.frame(chr="chr6",start=chromHMMTrack_43_tab[,2],
                                   end=chromHMMTrack_43_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_43_tab[,4],
                                   id=1:nrow(chromHMMTrack_43_tab),
                                   group=1:nrow(chromHMMTrack_43_tab))
chromHMMTrack_43_gr <- makeGRangesFromDataFrame(chromHMMTrack_43_mat, TRUE)
chromHMMTrack_43 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_43_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT43 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_43) <- colorchromHMM
displayPars(chromHMMTrack_43) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 

# data 44
chromHMMTrack_44_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_44_19_segments_AITD.bed")
chromHMMTrack_44_tab <- read.table(file=chromHMMTrack_44_path,header=FALSE,sep="\t")

chromHMMTrack_44_mat <- data.frame(chr="chr6",start=chromHMMTrack_44_tab[,2],
                                   end=chromHMMTrack_44_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_44_tab[,4],
                                   id=1:nrow(chromHMMTrack_44_tab),
                                   group=1:nrow(chromHMMTrack_44_tab))
chromHMMTrack_44_gr <- makeGRangesFromDataFrame(chromHMMTrack_44_mat, TRUE)
chromHMMTrack_44 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_44_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT44 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_44) <- colorchromHMM
displayPars(chromHMMTrack_44) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 

# data 45
chromHMMTrack_45_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_45_19_segments_AITD.bed")
chromHMMTrack_45_tab <- read.table(file=chromHMMTrack_45_path,header=FALSE,sep="\t")

chromHMMTrack_45_mat <- data.frame(chr="chr6",start=chromHMMTrack_45_tab[,2],
                                   end=chromHMMTrack_45_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_45_tab[,4],
                                   id=1:nrow(chromHMMTrack_45_tab),
                                   group=1:nrow(chromHMMTrack_45_tab))
chromHMMTrack_45_gr <- makeGRangesFromDataFrame(chromHMMTrack_45_mat, TRUE)
chromHMMTrack_45 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_45_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT45 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_45) <- colorchromHMM
displayPars(chromHMMTrack_45) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 

# data 86
chromHMMTrack_86_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_86_19_segments_AITD.bed")
chromHMMTrack_86_tab <- read.table(file=chromHMMTrack_86_path,header=FALSE,sep="\t")

chromHMMTrack_86_mat <- data.frame(chr="chr6",start=chromHMMTrack_86_tab[,2],
                                   end=chromHMMTrack_86_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_86_tab[,4],
                                   id=1:nrow(chromHMMTrack_86_tab),
                                   group=1:nrow(chromHMMTrack_86_tab))
chromHMMTrack_86_gr <- makeGRangesFromDataFrame(chromHMMTrack_86_mat, TRUE)
chromHMMTrack_86 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_86_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT86 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_86) <- colorchromHMM
displayPars(chromHMMTrack_86) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 

# data 87
chromHMMTrack_87_path<-paste0(dirFiles,"thyroid_chromHMM/CEMT_87_19_segments_AITD.bed")
chromHMMTrack_87_tab <- read.table(file=chromHMMTrack_87_path,header=FALSE,sep="\t")

chromHMMTrack_87_mat <- data.frame(chr="chr6",start=chromHMMTrack_87_tab[,2],
                                   end=chromHMMTrack_87_tab[,3],strand="*",
                                   score=0,feature=chromHMMTrack_87_tab[,4],
                                   id=1:nrow(chromHMMTrack_87_tab),
                                   group=1:nrow(chromHMMTrack_87_tab))
chromHMMTrack_87_gr <- makeGRangesFromDataFrame(chromHMMTrack_87_mat, TRUE)
chromHMMTrack_87 <- AnnotationTrack(genome="hg19",range=chromHMMTrack_87_gr,
                                    chromosome="chr6",
                                    name = "chromHMM CEMT87 cell-line",
                                    stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
displayPars(chromHMMTrack_87) <- colorchromHMM
displayPars(chromHMMTrack_87) <- list(rotation.title = 360, cex.title=0.5, 
                                      shape="box",stackHeight=0.8, 
                                      background.title="transparent", col.title="black", 
                                      col=NULL) 

listgviz_all <- c(list(itrack,gtrack, genetrack_highlighted,IgGsnp,
                   Immunesnp,AITDsnp),
              chromhmmPattern, 
              list(eQTL_blood),chromHMMTrack_40,chromHMMTrack_41,
              chromHMMTrack_42, chromHMMTrack_43,chromHMMTrack_44,
              chromHMMTrack_45,chromHMMTrack_86, chromHMMTrack_87,
              list(eQTL_thyroid))

displayPars(genetrack_highlighted) <- list(fontsize=14)
listgviz <- c(list(itrack,gtrack, genetrack_highlighted,IgGsnp,
                   Immunesnp,AITDsnp),
              chromhmmPattern, 
              list(eQTL_blood),chromHMMTrack_40,list(eQTL_thyroid))

## Check annotation tracks
plotTracks(listgviz, from=start, to=end, 
           fontsize=12, showId = TRUE,panel.only=TRUE)
dev.off()





##Independant immune cell traits based on 51 immune cell traits
source("multTest.r")
(meff_Li <- estimate.meff(immune_features_sub51_NA, method ="Li"))
(meff_Cheverud <- estimate.meff(immune_features_sub51_NA, method ="Cheverud"))

## select only 51 immune cell traits for Association analysis between Immune cell traits and AITD or TPOAb level
AITDimmunefile <-"ImmuneTrait_AITD_20160221_all_results_withlist_immuneTraitssignificantIgG.csv"
AITDimmune <-read.table(file=AITDimmunefile,header=TRUE,sep=",")
dim(AITDimmune)
AITDimmune51 <- AITDimmune[which(AITDimmune$immuneTrait %in% rownames(immune_features_sub51)),]
dim(AITDimmune51)
AITDimmunefilesub51 <-"ImmuneTrait_AITD_20160221_all_results_withlist_immuneTraitssignificantIgG_sub51.csv"
write.csv(AITDimmune51,file=AITDimmunefilesub51,row.names=FALSE)

TPOAbimmunefile <-"ImmuneTrait_antiTPO_20160221_all__results_final.csv"
TPOAbimmune <-read.table(file=TPOAbimmunefile,header=TRUE,sep=",")
dim(TPOAbimmune)
TPOAbimmune51 <- TPOAbimmune[which(TPOAbimmune$immuneTrait %in% rownames(immune_features_sub51)),]
dim(TPOAbimmune51)
TPOAbimmunefilesub51 <-"ImmuneTrait_antiTPO_20160221_all__results_immuneTraitssignificantIgG_sub51.csv"
write.csv(TPOAbimmune51,file=TPOAbimmunefilesub51,row.names=FALSE)

##### to fill up the supplementarty tables
###
#General info on immuneWAS- glycome
immuneGlufile <- "/Volumes/DATABAR/immune_traits/glycan_proteins_metadata_residual.txt"
immuneGlu <-read.table(file=immuneGlufile,header=T,sep="\t")
dim(immuneGlu)

mean(immuneGlu$Age)
sd(immuneGlu$Age)

## 
# General info on immuneWAS-AITD
name_fileRData <- "ds.age.and.batch.corrected.quantile.normalised_all.RData"
load(name_fileRData)

immune_info <- read.table(file="immuno_poppante_AGE_BATCH.covar",sep="\t")
colnames(immune_info)[2] <- "ID"
colnames(immune_info)[3] <- "age_immune"

data0 <-merge(immune_info[,c(2,3)],ds.age.and.batch.corrected.quantile.normalised_interesting_date_clean,by="ID")

AITD <- read.table(file="AITD_Roche_Abbott.txt",sep="\t",header=T)

data1_tmp <-merge(AITD,data0,by="ID")
dim(data1_tmp)

library(plyr)
(max_le_by_cont <-
    ddply(data1_tmp, ~ batch + status_AITD , summarize, num=length(status_AITD),
          meanAge = mean(age_ROCHE),sdAge =sd(age_ROCHE)))


roche <- read.table(file='roche/TFT_materList_FIN_Roche2011_GL_allphenotype_final_AITD.txt',sep="\t",header=T)
abbot <- read.table(file='abbott/Full_Data_Set_Thyroid_abbot2011_GL_allphenotype_final_AITD.txt',sep="\t",header=T)


data1_tmpa <- data1_tmp[which(data1_tmp$batch %in% "Abbot"),]
data1_tmpr <- data1_tmp[which(data1_tmp$batch %in% "roche"),]


data1_tmpr1 <- merge(data1_tmpr,roche,by="ID")
dim(data1_tmpr1)

(max_le_by_cont <-
    ddply(data1_tmpr1 ~ status_AITD.y , summarize, num=length(status_AITD.x),
          meanAge = mean(age_ROCHE.y),sdAge =sd(age_ROCHE.y),
          meanTPOAb = mean(AntiTPO_Roche),sdTPOAb =sd(AntiTPO_Roche)))
mean(data1_tmpr1[which(data1_tmpr1$status_AITD.x == 1),"AntiTPO_Roche"])
sd(data1_tmpr1[which(data1_tmpr1$status_AITD.x == 1),"AntiTPO_Roche"])
mean(data1_tmpr1[which(data1_tmpr1$status_AITD.x == 2),"AntiTPO_Roche"])
sd(data1_tmpr1[which(data1_tmpr1$status_AITD.x == 2),"AntiTPO_Roche"])

data1_tmpa1 <- merge(data1_tmpr,abbot,by="ID")
mean(data1_tmpa1[which(data1_tmpa1$status_AITD.x == 1),"AntiTPO_abbot"])
sd(data1_tmpa1[which(data1_tmpa1$status_AITD.x == 1),"AntiTPO_abbot"])
mean(data1_tmpa1[which(data1_tmpa1$status_AITD.x == 2),"AntiTPO_abbot"])
sd(data1_tmpa1[which(data1_tmpa1$status_AITD.x == 2),"AntiTPO_abbot"])
