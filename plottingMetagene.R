setwd("~/Documents/sequencingData/20170323_dSMF_N2/dSMFseq_scripts")
source('./callAllCs.r') #to call all Cs
source('./useful_functionsV1.r')
source("~/Documents/MeisterLab/dSMF/TSS/scripts/getTSS_functions.R")
source("~/Documents/MeisterLab/dSMF/usefulFunctions.R")
source("~/Documents/MeisterLab/dSMF/usefulFunctions.R")
source("~/Documents/MeisterLab/dSMF/findMotifs/findMotifs_functions.R")
library(Biostrings)
library(Gviz)

#metageneMatrix created by dSMFseqAnalysis_3plot1.R


designT=read.csv('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/finalChosenList_tss_dc_tissues.csv',
                 stringsAsFactors=FALSE,header=TRUE)
design=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/amplicons_stranded.bed')
tss=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/ampliconTSS.bed')

amplicons=design
amplicons=resize(amplicons,500,fix='center')

TSS=tss
TSS2=resize(tss,500,fix='center')


methMatrices_tss<-readRDS("./allMethMatrices_tss.RDS")
metaGene<-matrix(data=NaN,nrow=length(methMatrices_tss),ncol=500)
for (i in sl(methMatrices_tss)) {
  st=TSS[i]
  mat<-methMatrices_tss[[i]]
  if (!is.null(dim(mat))) {
    if(dim(mat)[1]>20) {
      pos<-as.numeric(colnames(mat))
      if(as.character(strand(st))=="+") {
        startWin<-start(st)-250
        startPrhc<-pos-startWin
      } else {
        startWin<-start(st)+250
        startPrhc<-startWin-pos
      }
      frMeth<-colMeans(mat,na.rm=TRUE)
      #matRow<-rep(NaN,500)
      #matRow[startPrhc]<-frMeth
      metaGene[i,startPrhc]<-frMeth
    }
  }
}

pdf("metagenePlot.pdf",paper="a4",height=11,width=8)
par(mfrow=c(2,1))
plot(c(1:500),colMeans(metaGene,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="fraction methylated",
     main="Average fraction methylated around TSS")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess(colMeans(metaGene,na.rm=TRUE)~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")

plot(c(1:500),1-colMeans(metaGene,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="fraction protected",
     main="Average fraction protected around TSS")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess((1-colMeans(metaGene,na.rm=TRUE))~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")
dev.off()


#simpsons paradox is a big deal with methylation profiles because many positions might not
#have a CG to methylate, so a  0 there is not the same as 0 methylation. see:
#https://liorpachter.wordpress.com/2015/07/13/how-to-average-genome-wide-data/
#write metaGene file for perl COMPARE scipt
#remove rows with no data
mg<-metaGene[-which(rowSums(is.na(metaGene))==500),]
#write as tab separated file
write.table(mg,file="metaGeneMat.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,na="NaN")
# run COMPARE.pl on command line. read in the imputed matrix:
mg1  <-read.table("./imputedMetagene/impMG/Inference_Matrix.all_inf")


pdf("imp_metagenePlot.pdf",paper="a4",height=11,width=8)
par(mfrow=c(2,1))
plot(c(1:500),colMeans(mg1,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="fraction methylated (imputed)",
     main="Average fraction methylated around TSS")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess(colMeans(mg1,na.rm=TRUE)~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")

plot(c(1:500),1-colMeans(mg1,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="fraction protected (imputed)",
     main="Average fraction protected around TSS")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess((1-colMeans(mg1,na.rm=TRUE))~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")
dev.off()
#plots look roughly similar as simple metagenes, but perhaps a bit flattened, and perhaps a bit
#more of a bump at -50?


#importing gff3 data
library("rtracklayer")
fileName<-"~/Documents/MeisterLab/Datasets/nucleosome_GuFire2010_GSM432720/GSM432720_totalNucleosome_N2_seq_WS235.bed"
#fileName<-list.files(path)[1]
GuFireNuc<-import.bed(fileName)
Nuc<-subsetByOverlaps(GuFireNuc, TSS2)

nucOccup<-matrix(data=0,nrow=96,ncol=500)
for (st in sl(TSS2)) {
  if (as.character(strand(TSS2[st]))=="+") {
    matRow<-rep(0,500)
    for(i in 1:500) {
      tempGR<-GRanges(seqnames=seqnames(TSS2[st]),IRanges(start=start(TSS2[st]+i-1),width=1),strand=strand(TSS2[st]))
      matRow[i]<-countOverlaps(tempGR,Nuc)
    }
  } else {
    matRow<-rep(0,500)
    for(i in 1:500) {
      tempGR<-GRanges(seqnames=seqnames(TSS2[st]),IRanges(start=start(TSS2[st]+i-1),width=1),strand=strand(TSS2[st]))
      matRow[i]<-countOverlaps(tempGR,Nuc)
    }
    matRow<-rev(matRow)
  }
  nucOccup[st,]<-matRow
}

plot(c(1:500),colMeans(nucOccup,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="nucleosome occupancy",
     main="Average nucleosome occupancy around TSS (Gu et al. 2010)",cex.main=0.9)
#barplot(nucOccup,col="grey")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess((colMeans(nucOccup,na.rm=TRUE))~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")

#### L3 data in gff3 format
path="~/Documents/MeisterLab/Datasets/modEncode_44_174_599_3276_2764_2763_3564_2407_2408_2409_3206_3562/modEncode_3562/interpreted_data_files/"
fileName<-"H3_L3_ChIP-seq_R1_WS220_3562_catted_1.gff3"
import.gff3(paste0(path,fileName))
#error saying too many NAs, even when comment lines removed


######################### Adults MNase data
#create seqinfo object for conversion of wig to bigwig
library(BSgenome.Celegans.UCSC.ce10)
SI<-seqinfo(BSgenome.Celegans.UCSC.ce10)
#SI<-ucsc2wb(SI)
genome(SI)<-"ce10" #"WS220"
#convert wig to bigwig
path="~/Documents/MeisterLab/Datasets/MNase_Ad_Steiner2012_GSM777719/"
fileName<-"GSM777719_Total_nuclei_20110313_3.UCSC.wig"
#wigToBigWig(paste0(path,fileName),SI)
AdNuc<-import(paste0(path,"GSM777719_Total_nuclei_20110313_3.UCSC.bw"),format="BigWig",
                which=TSS2)
nucOcc<-matrix(data=0,nrow=96,ncol=500)
for (st in sl(TSS2)) {
  currTSS<-TSS2[st]
  strnd<-as.character(strand(currTSS))
  tempGR<-subsetByOverlaps(AdNuc,TSS2[st])
  matRow=rep(0,500)
  if (strnd=="+") {
    allStarts<-start(tempGR)-start(currTSS)+1
    allEnds<-end(tempGR)-start(currTSS)+1
    for (i in sl(tempGR)) {
      matRow[allStarts[i]:allEnds[i]]<-mcols(tempGR[i])$score
    }
  } else {
    allStarts<-start(tempGR)-start(currTSS)+1
    allEnds<-end(tempGR)-start(currTSS)+1
    for (i in sl(tempGR)) {
      matRow[allStarts[i]:allEnds[i]]<-mcols(tempGR[i])$score
    }
    matRow<-rev(matRow)
  }
  nucOcc[st,]<-matRow
}

pdf("metaGenePlotsWnucleosomes.pdf",paper="a4",height=11,width=8)

par(mfrow=c(2,1))
plot(c(1:500),1-colMeans(metaGene,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="fraction protected",
     main="Average fraction of protected cytosines around TSS (dSMF)",cex.main=1.2)
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess((1-colMeans(metaGene,na.rm=TRUE))~c(1:500),span=0.3)
lines(predict(l),lwd=3)
abline(v=251,col="red")

plot(c(1:500),colMeans(nucOcc,na.rm=TRUE),type="l",lwd=3,col="black",xaxt="n",
     xlab="position relative to TSS",ylab="nucleosome occupancy",
     main="Average nucleosome occupancy around TSS (Steiner et al. 2012)",cex.main=1.2)
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
lines(predict(l))
abline(v=251,col="red")

dev.off()


######################### Embryos MNase data
#create seqinfo object for conversion of wig to bigwig
library(BSgenome.Celegans.UCSC.ce10)
SI<-seqinfo(BSgenome.Celegans.UCSC.ce10)
#SI<-ucsc2wb(SI)
genome(SI)<-"ce10" #"WS220"
#convert wig to bigwig
path="~/Documents/MeisterLab/Datasets/MNase_emb_Ooi2010_GSM514735/"
fileName<-"GSM514735_Embryos_Rep2_RawCoverage.wig"
wigToBigWig(paste0(path,fileName),SI)
embNuc<-import(paste0(path,"GSM514735_Embryos_Rep2_RawCoverage.bw"),format="BigWig",
              which=TSS2)
nucOccEmb<-matrix(data=0,nrow=96,ncol=500)
for (st in sl(TSS2)) {
  currTSS<-TSS2[st]
  strnd<-as.character(strand(currTSS))
  tempGR<-subsetByOverlaps(embNuc,TSS2[st])
  matRow=rep(0,500)
  if (strnd=="+") {
    allStarts<-start(tempGR)-start(currTSS)+1
    allEnds<-end(tempGR)-start(currTSS)+1
    for (i in sl(tempGR)) {
      matRow[allStarts[i]:allEnds[i]]<-mcols(tempGR[i])$score
    }
  } else {
    allStarts<-start(tempGR)-start(currTSS)+1
    allEnds<-end(tempGR)-start(currTSS)+1
    for (i in sl(tempGR)) {
      matRow[allStarts[i]:allEnds[i]]<-mcols(tempGR[i])$score
    }
    matRow<-rev(matRow)
  }
  nucOccEmb[st,]<-matRow
}



pdf("metaGenePlotsWnucleosomes.pdf",paper="a4",height=11,width=8)

par(mfrow=c(3,1))
plot(c(1:500),1-colMeans(metaGene,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="fraction protected",
     main="Average fraction protected around TSS")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess((1-colMeans(metaGene,na.rm=TRUE))~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")

plot(c(1:500),colMeans(nucOcc,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="nucleosome occupancy",
     main="Average nucleosome occupancy around TSS (Steiner et al. 2012)",cex.main=0.9)
#barplot(nucOccup,col="grey")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess((colMeans(nucOccup,na.rm=TRUE))~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")

plot(c(1:500),colMeans(nucOccEmb,na.rm=TRUE),pch=16,cex=0.8,col="#77777755",xaxt="n",
     xlab="position relative to TSS",ylab="nucleosome occupancy",
     main="Average nucleosome occupancy around TSS (Ooi et al. 2010)",cex.main=0.9)
#barplot(nucOccup,col="grey")
axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
l<- loess((colMeans(nucOccup,na.rm=TRUE))~c(1:500),span=0.3)
lines(predict(l))
abline(v=251,col="red")

dev.off()



#TSSregions<-paste0(seqnames(TSS2),":",start(TSS2),"-",end(TSS2))
#write.table(TSSregions,"TSSregions.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
