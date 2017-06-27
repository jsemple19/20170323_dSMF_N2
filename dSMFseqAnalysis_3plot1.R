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
amplicons=resize(amplicons,500,fix='center')

TSS=tss
TSS2=resize(tss,500,fix='center')
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

#amplicon<-amplicons[84]

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
pdf("./methMatPlotsAK.pdf",paper="a4",height=11,width=8)
par(mfrow=c(2,1))
for (i in sl(TSS)) {
   st=TSS[i]
   mat<-tryCatch(
      { mat<-getMethMatrix(amplicons,designT,i,NOMEproj,"N2_DE_ampV001") },
      error=function(cond) {return(NULL)})
   if (length(mat)!=0) {
      vRhc=VectorizeReads(st, mat)
   #plot the footprints
      BR=c('black','#E6E6E6')
      colors=BR
      if(as.character(strand(st))=="+") {
        startPrhc<-vRhc[[1]]-start(st)
      } else {
        startPrhc<-start(st)-vRhc[[1]]
      }
      plot(rev(startPrhc),rev(vRhc[[2]]),pch='_',cex=1,col=colors[as.factor(rev(vRhc[[3]]))],xlim=c(-150,150),
           xlab="position relative to TSS",ylab="molecules",
           main=paste0(designT[i,]$Gene_WB_ID,", expn decile: ",designT[i,]$decile,", tissues expd: ",
                       designT[i,]$SpencerLarval11,", ",strand(st)),cex.main=0.7)
      legend("bottomright",legend=c("protected","open (meth)"),col=colors,lty=1,lwd=2,cex=0.6,
             box.col="white")
   }
}
dev.off()



#mat=extractMultipleMatrices(resize(st,300,fix='center'),AMPproj,
#                            as.character(spNames.amp[sp.sbs[5]]),Dmelanogaster,1)
#vectorise the read information with
getMethMatrix<-function(amplicons,designT,ampNum,proj,sampleName,genome=Celegans) {
   amplicon<-amplicons[ampNum]
   geneName<-mcols(amplicon)$name
   chr<-strsplit(designT[designT$Gene_WB_ID==geneName,]$Amplicon,":")[[1]][1]
   tss<-designT[designT$Gene_WB_ID==geneName,]$maxTSS
   allCs<-getCMethMatrix(proj,amplicon,sampleName)
   onlyCG_GC<-getGCMatrix(allCs,chr=chr,genome=genome,conv.rate=80,destrand=FALSE)
   mCG_GC<-mergeGC_CGmat(onlyCG_GC)
   #reads<-dim(mCG_GC)[1]
   mCG_GC.o<-mCG_GC[hclust(dist(mCG_GC))$order,]
   return(mCG_GC.o)
}

VectorizeReads=function(target,mergedMat){ #vectorize reads from matrices in order to plot in base space

  # 			 mergedMat=matList[[i]][naidx,][h$rowInd,]
  GCmet=as.vector(t(mergedMat))
  GCstartP=rep(as.numeric(colnames(mergedMat)),nrow(mergedMat))
  GCreadID=unlist(lapply(seq(nrow(mergedMat)),function(i){rep(i,ncol(mergedMat))}))
  GCgr=GRanges(seqnames(target),IRanges(as.numeric(colnames(mergedMat)),as.numeric(colnames(mergedMat))))
  list(GCstartP,GCreadID,GCmet)
}



extractAllMethMat<-function(TSS,amplicons,designT,proj,sampleName) {
  allMethMat<-list()
  for (i in sl(TSS)) {
    st=TSS[i]
    mat<-tryCatch(
      { mat<-getMethMatrix(amplicons,designT,i,proj,sampleName) },
      error=function(cond) {return(NULL)})
    if (is.null(mat)) {
      allMethMat[[designT[i,]$Gene_WB_ID]]<-"0 reads"
    } else {
      allMethMat[[designT[i,]$Gene_WB_ID]]<-mat
    }
  }
  return(allMethMat)
}

methMatrices<-extractAllMethMat(TSS,amplicons,designT,NOMEproj,"N2_DE_ampV001")
saveRDS(methMatrices,"./allMethMatrices.RDS")

unlist(lapply(lapply(methMatrices_tss,dim), '[[',1))

methMatrices_tss<-readRDS("./allMethMatrices_tss.RDS")
###############################################################
pdf("./methMatPlotsAKgt200.pdf",paper="a4",height=11,width=8)
par(mfrow=c(4,2))
plotMethMatrices(TSS,amplicons,designT,methMatrices_tss)
dev.off()

plotMethMatrices<-function(TSS,amplicons,designT,methMatrices) {
  for (i in sl(methMatrices)) {
    st=TSS[i]
    mat<-methMatrices[[i]]
    if (!is.null(dim(mat))) {
        if(dim(mat)[1]>200) {
          vRhc=VectorizeReads(st, mat)
      #plot the footprints
          BR=c('black','#D5D5D5')
          colors=BR
          if(as.character(strand(st))=="+") {
            startPrhc<-vRhc[[1]]-start(st)
        } else {
          startPrhc<-start(st)-vRhc[[1]]
        }
        #molecules<-dim(mat)[1]
      #oldmai<-par()$mai
        #newmai<-oldmai+c(2-2*molecules/600,0,0,0)
        #par(mai=newmai)
        plot(rev(startPrhc),rev(vRhc[[2]]),pch='_',cex=1,col=colors[as.factor(rev(vRhc[[3]]))],xlim=c(-250,250),
           xlab="position relative to TSS",ylab="molecules",
           main=paste0(designT[i,]$Gene_WB_ID,", expn decile: ",designT[i,]$decile,", tissues expd: ",
                       designT[i,]$SpencerLarval11,", ",strand(st)),cex.main=0.7)
        abline(v=0,col="red")
        legend("bottomright",legend=c("protected","accessible"),col=colors,pch="_",lwd=2,cex=0.7,
            bg="white")
      #par(mai=oldmai)
      }
    }
  }
}

############# create metagene matrix +-250 around TSS
TSS2=resize(tss,500,fix='center')
mcols(TSS2)$name<-mcols(amplicons)$name
methMatrices_tss<-extractAllMethMat(TSS,TSS2,designT,NOMEproj,"N2_DE_ampV001")
methMatrices_tss<-readRDS("./allMethMatrices_tss.RDS")


############# plot metagene
methMatrices_tss<-readRDS("./allMethMatrices_tss.RDS")
metaGene<-matrix(data=NaN,nrow=length(methMatrices_tss),ncol=500)
row.names(metaGene)<-names(methMatrices_tss)
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
#row.names(metaGene)<-names(methMatrices_tss)
#mg<-metaGene[-which(rowSums(is.na(metaGene))==500),]

pdf("fractionProtected_line.pdf",paper="a4",height=11,width=8)
par(mfrow=c(4,2))
for (i in 1:dim(mg)[1]) {
  geneName=row.names(mg)[i]
  Ti=which(designT$Gene_WB_ID==geneName) #correct indexing to access data table
  plot(c(1:500),1-mg[i,],pch=16,xaxt="n", ylim=c(0,1),
       xlab="position relative to TSS",ylab="fraction protected",
       main=paste0(designT[Ti,]$Gene_WB_ID,", expn decile: ",designT[Ti,]$decile,", tissues expd: ",
              designT[Ti,]$SpencerLarval11,", ",designT[Ti,]$strand),cex.main=0.7)
  j<-!is.na(mg[i,])
  lines(c(1:500)[j],1-mg[i,j])
  axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
  abline(v=251,col="red")
}
dev.off()

########### create cartoons of motifs

motifs<-readRDS("~/Documents/MeisterLab/dSMF/findMotifs/amplicons_motifs.RDS")
pdf("motifCartoons.pdf",paper="a4",height=11,width=8)
par(mfrow=c(4,2))
for (i in sl(motifs)) {
  plot(c(1,500),c(0.5,0.5),ylim=c(0,1),xaxt="n",yaxt="n",type="l",ylab="",bty="n",
       xlab="position relative to TSS", lwd=3, cex.main=0.7,
       main=paste0(designT[i,]$Gene_WB_ID,", expn decile: ",designT[i,]$decile,", tissues expd: ",
                  designT[i,]$SpencerLarval11,", ",designT[i,]$strand))
  lines(c(250,250),c(0.5,0.65),col="black",lwd=2)
  arrows(250,0.65,290,0.65,lwd=2,code=2,length=0.13)
  inr<-as.numeric(strsplit(mcols(motifs[i,])$inr_relTSSmax,"[:-]")[[1]][2:3])
  tata<-as.numeric(strsplit(mcols(motifs[i,])$tata_relTSSmax,"[:-]")[[1]][2:3])
  cgrich<-as.numeric(strsplit(mcols(motifs[i,])$cgrich_relTSSmax,"[:-]")[[1]][2:3])
  if(as.character(mcols(motifs[i,])$strand)=="+") {
    inr<-inr-mcols(motifs[i,])$maxTSS
    tata<-tata-mcols(motifs[i,])$maxTSS
    cgrich<-cgrich-mcols(motifs[i,])$maxTSS
  } else {
    inr<-mcols(motifs[i,])$maxTSS-inr
    tata<-mcols(motifs[i,])$maxTSS-tata
    cgrich<-mcols(motifs[i,])$maxTSS-cgrich

  }
  polygon(c(inr,rev(inr))+250,c(0.43,0.43,0.57,0.57),col="red")
  text(inr[1]+251,0.35,labels="Inr",col="red",adj=0.1)
  polygon(c(tata,rev(tata))+250,c(0.43,0.43,0.57,0.57),col="blue")
  text(tata[1]+251,0.35,labels="TATA",col="blue",adj=0.65)
  polygon(c(cgrich,rev(cgrich))+250,c(0.43,0.43,0.57,0.57),col="orange")
  text(cgrich[1]+251,0.35,labels="CG rich",col="orange",adj=1)
  axis(1,at=seq(1,501,50),labels=seq(-250,250,50))
}
dev.off()



hc=hclust(dist(s.mat))
vRhc=VectorizeReads(st, s.mat[hc$order,])

vRhc=VectorizeReads(st, a)

#plot the footprints
BR=c('black','#E6E6E6')
colors=BR

startPrhc<-vRhc[[1]]-start(st)

plot(rev(startPrhc),rev(vRhc[[2]]),pch='_',cex=1,col=colors[as.factor(rev(vRhc[[3]]))],xlim=c(-150,150))

plot(rev(vRhc[[1]]),rev(vRhc[[2]]),pch='_',cex=1,col=colors[as.factor(rev(vRhc[[3]]))],xlim=c(-150,150))
