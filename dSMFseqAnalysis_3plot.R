library(QuasR)
library("BSgenome.Celegans.UCSC.ce11")

setwd("~/Documents/sequencingData/20170323_dSMF_N2/dSMFseq_scripts")
source('./callAllCs.r') #to call all Cs
source('./useful_functionsV1.r') #load the ranges



path='~/Documents/sequencingData/20170323_dSMF_N2/'

NOMEproj=qAlign(sampleFile='./sample_BAM.tmp', 
                genome="BSgenome.Celegans.UCSC.ce11", 
                paired="fr", 
                bisulfite="dir", 
                projectName="dSMF_N2", 
                clObj=cluObj)
NOMEaln=as.data.frame(alignments(NOMEproj)[[1]])

samples=NOMEaln$SampleName



#PolII_in<-readRDS(paste0(path,"/NomeAmpliconCEmeth.rds"))

designT=read.csv('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/finalChosenList_tss_dc_tissues.csv',
                 stringsAsFactors=FALSE,header=TRUE)
design=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/amplicons_stranded.bed')
tss=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/ampliconTSS.bed')

#maxTSS=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/maxTSS.bed')
#start(maxTSS)<-end(maxTSS)

amplicons=design
amplicons=resize(amplicons,600,fix='center')

TSS=tss
TSS2=resize(TSS,500,fix='center')
strand(TSS2)='+'

#meth_gr<-readRDS(paste0(path,'/Amplicon_raw_methylation_PolII_covered_amplicons.rds'))

#cO=10
#PolII_in=call_context_methylation(meth_gr,cO,Celegans)

#PolII_inm=unlist(GRangesList(PolII_in))
# saveRDS(PolII_inm,'/work2/gschub/arnaud/nomeseq/WG/DM/amplicon/methylation_calls/PolII_inhib_av10.rds')
#saveRDS(PolII_inm,paste0(path,'PolII_inhib_av10.rds'))
#naix=apply(as.matrix(elementMetadata(PolII_inm[,1])),1,function(x){sum(is.na(x))==1})
#PolII_inm<-PolII_inm[!naix]
#saveRDS(PolII_inm,paste0(path,'PolII_inhib_av10.rds'))

PolII_inm<-readRDS(paste0(path,'PolII_inhib_av10.rds'))

amplicon<-amplicons[84]

plotMethPos<-function(amplicon,NOMEproj,designT,genome=Celegans) {
   geneName<-mcols(amplicon)$name
   chr<-strsplit(designT[designT$Gene_WB_ID==geneName,]$Amplicon,":")[[1]][1]
   TSS<-designT[designT$Gene_WB_ID==geneName,]$maxTSS
   allCs<-getCMethMatrix(NOMEproj,amplicon,"N2_DE_ampV001")
   onlyCG_GC<-getGCMatrix(allCs,chr=chr,genome=Celegans,conv.rate=80,destrand=FALSE)
   mCG_GC<-mergeGC_CGmat(onlyCG_GC)
   methPos<-as.integer(colnames(mCG_GC))
   mCG_GC.d<-dist(mCG_GC,method = "euclidian",upper=TRUE)
   mCG_GC.hc<-hclust(mCG_GC.d)
   mCG_GC.hco<-mCG_GC[mCG_GC.hc$order,]
   image(t(mCG_GC.hco),axes=FALSE,col=grey(seq(0,0.8,length=2)),ylab="molecule",
      xlab="methylation site") #plots 0 as black and 1 as grey
   i<-max(which(methPos<TSS))
   abline(v=i/length(methPos),col="red",lwd=2)
   #text(i/length(methPos),1.05,labels="TSS",col="red")
   title(geneName)
}
#remove<-seq(7,96,by=12)
#goodAmps<-amplicons[-remove]
pdf("plotMethPos.pdf",paper="a4",height=11,width=8)
par(mfrow=c(2,2))
for (a in 1:length(amplicons)) {
   try(plotMethPos(amplicons[a],NOMEproj,designT),silent=TRUE)
   
}
dev.off()


###########
library(Gviz)
library(GenomicFeatures)
txdb<-makeTxDbFromGFF("~/Documents/MeisterLab/GenomeVer/annotations/c_elegans.PRJNA13758.WS250.annotations.gff3")
genes<-genes(txdb)
gtrack<-GenomeAxisTrack()
atrack<-AnnotationTrack(PolII_inm,name="methylation")
#grtrack<-GeneRegionTrack(genes,)

plotTracks(list(gtrack,atrack))
amplicon<-amplicons[84]

ol<-findOverlaps(amplicon,PolII_inm,ignore.strand=TRUE)

thisAmp<-PolII_inm[subjectHits(ol)]
thisAmp<-sort(thisAmp)
df<-data.frame(thisAmp)
df$NAcols<-"ok"
df$NAcols[naCols]<-"NA"
write.table(df,"WBGene00016684sites.txt")

pdf("fractionProtected.pdf",paper="a4",height=11,width=8)
par(mfrow=c(2,1))
plot(start(thisAmp),1-mcols(thisAmp)$V1,type="p",pch=16,xlim=c(start(amplicon),end(amplicon)),
     col=as.factor(mcols(thisAmp)$type),xlab="genomic position",ylab="fraction protected")
lines(start(thisAmp),1-mcols(thisAmp)$V1)
legend("topleft",legend=levels(as.factor(mcols(thisAmp)$type)),fill=c(1,2,3))
title(mcols(amplicon)$name)
abline(v=start(maxTSS[mcols(maxTSS)$name=="WBGene00016684"]),col="red",lwd=2)

plot(start(thisAmp),1-mcols(thisAmp)$V1,type="p",pch=16,xlim=c(start(amplicon),end(amplicon)),
     xlab="genomic position",ylab="fraction protected")
lines(start(thisAmp),1-mcols(thisAmp)$V1)
title(mcols(amplicon)$name)
abline(v=start(maxTSS[mcols(maxTSS)$name=="WBGene00016684"]),col="red",lwd=2)
dev.off()


pdf("fractionProtected_all.pdf",paper="a4",height=11,width=8)
par(mfrow=c(3,2))
for (i in 1:length(amplicons)){
   amplicon<-amplicons[i]
   plotFractionProtected(amplicon,PolII_inm,designT)
}
dev.off()
genome=Celegans
plotFractionProtected<-function(amplicon,C_GRanges,designT) {
   geneName<-mcols(amplicon)$name
   chr<-strsplit(designT[designT$Gene_WB_ID==geneName,]$Amplicon,":")[[1]][1]
   maxTSS<-designT[designT$Gene_WB_ID==geneName,]$maxTSS
   strand<-designT[designT$Gene_WB_ID==geneName,]$strand
   decile<-designT[designT$Gene_WB_ID==geneName,]$decile
   ol<-findOverlaps(amplicon,PolII_inm,ignore.strand=TRUE)
   if (length(ol)==0){
      return(1)
   }
   thisAmp<-PolII_inm[subjectHits(ol)]
   thisAmp<-thisAmp[order(start(thisAmp))]
   plot(start(thisAmp),1-mcols(thisAmp)$V1,type="p",pch=16,xlim=c(start(amplicon),end(amplicon)),
     col=as.factor(mcols(thisAmp)$type),xlab="genomic position",ylab="fraction protected")
   lines(start(thisAmp),1-mcols(thisAmp)$V1)
   legend("topleft",legend=levels(as.factor(mcols(thisAmp)$type)),fill=c(1,2,3))
   title(paste(geneName, "strand:",strand,"expnDecile:",decile))
   abline(v=maxTSS,col="blue",lwd=2)
}


amplicon<-amplicons[84]
geneName<-mcols(amplicon)$name
chr<-strsplit(designT[designT$Gene_WB_ID==geneName,]$Amplicon,":")[[1]][1]
maxTSS<-designT[designT$Gene_WB_ID==geneName,]$maxTSS
allCs<-getCMethMatrix(NOMEproj,amplicon,"N2")
onlyCG_GC<-getGCMatrix(allCs,chr=chr,genome=genome,conv.rate=80,destrand=FALSE)
mCG_GC<-mergeGC_CGmat(onlyCG_GC)
reads<-dim(mCG_GC)[1]
#naCols<-which(colSums(is.na(mCG_GC))==reads)
#mCG_GC<-mCG_GC[,-naCols]
#mCG_GC<-na.omit(mCG_GC)

dms<-c("euclidean","manhattan","binary","minkowski")
cms<-c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")

pdf("clusteringParameters.pdf",paper="a4",height=11,width=8)
for (dm in dms) {
   for (cm in cms) {
      methPos<-as.integer(colnames(mCG_GC))
      distanceMeasure<-dm
      mCG_GC.d<-dist(mCG_GC,method = distanceMeasure ,upper=FALSE)
      clustMethod<-cm
      mCG_GC.hc<-hclust(mCG_GC.d,method=clustMethod)
      mCG_GC.hco<-mCG_GC[mCG_GC.hc$order,]
      image(t(mCG_GC.hco),axes=FALSE,col=grey(seq(0,0.8,length=2)),ylab="molecule",
         xlab="methylation site")
      i<-max(which(methPos<maxTSS))
      abline(v=i/length(methPos),col="red",lwd=2)
      #text(i/length(methPos),1.05,labels="TSS",col="red")
      title(paste(geneName,distanceMeasure,clustMethod))
   }
}
dev.off()



amplicon<-amplicons[84]
geneName<-mcols(amplicon)$name
chr<-strsplit(designT[designT$Gene_WB_ID==geneName,]$Amplicon,":")[[1]][1]
maxTSS<-designT[designT$Gene_WB_ID==geneName,]$maxTSS
allCs<-getCMethMatrix(NOMEproj,amplicon,"N2")
onlyCG_GC<-getGCMatrix(allCs,chr=chr,genome=genome,conv.rate=80,destrand=FALSE)
mCG_GC<-mergeGC_CGmat(onlyCG_GC)
reads<-dim(mCG_GC)[1]
methPos<-as.integer(colnames(mCG_GC))
distanceMeasure<-"euclidian"
mCG_GC.d<-dist(mCG_GC,method = distanceMeasure ,upper=FALSE)
clustMethod<-"complete"
mCG_GC.hc<-hclust(mCG_GC.d,method=clustMethod)
mCG_GC.hco<-mCG_GC[mCG_GC.hc$order,]
plot(x=c(start(amplicon),end(amplicon)),y=c(0,reads),type="n",xlab="position",ylab="molecule")
for (rw in 1:reads) {
   for (cl in 1:dim(mCG_GC.hco)[2]) {
      lines(x=c(rw,rw+1),y=c(methPos[cl],methPos[cl]),col=mCG_GC.hco[rw,cl],lwd=2)
   } 
}
image(t(mCG_GC.hco),axes=FALSE,col=grey(seq(0,0.8,length=2)),ylab="molecule",
      xlab="methylation site")
i<-max(which(methPos<maxTSS))
abline(v=i/length(methPos),col="red",lwd=2)
#text(i/length(methPos),1.05,labels="TSS",col="red")
title(paste(geneName,distanceMeasure,clustMethod))


###################3 single molecule plotting from Arnaud
#you can extract data with
st= GR()#region of intest (i.e TSS GRanges)
for (i in TSS) {
   st=TSS[i]
   
}
st=TSS[84]
mat=extractMultipleMatrices(resize(st,300,fix='center'),AMPproj,
                            as.character(spNames.amp[sp.sbs[5]]),Dmelanogaster,1)
#vectorise the read information with

extractAllMethMat<-function()

hc=hclust(dist(s.mat))
vRhc=VectorizeReads(st, s.mat[hc$order,])

vRhc=VectorizeReads(st, a)

#plot the footprints
BR=c('black','#E6E6E6')
colors=BR

startPrhc<-vRhc[[1]]-start(st)

plot(rev(startPrhc),rev(vRhc[[2]]),pch='_',cex=1,col=colors[as.factor(rev(vRhc[[3]]))],xlim=c(-150,150))

plot(rev(vRhc[[1]]),rev(vRhc[[2]]),pch='_',cex=1,col=colors[as.factor(rev(vRhc[[3]]))],xlim=c(-150,150))
