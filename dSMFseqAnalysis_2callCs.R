#/work/gbioinfo/Appz/R/R-devel/bin/R

library(QuasR)
setwd("~/Documents/sequencingData/20170323_dSMF_N2/dSMFseq_scripts")

source('./callAllCs.r') #to call all Cs
source('./useful_functionsV1.r') #load the ranges

library("BSgenome.Celegans.UCSC.ce11")



#load states ranges

   #all genome

   #process the samples in chunks

   path='~/Documents/sequencingData/20170323_dSMF_N2/'
   #my.alignmentsDir=paste(path,'aln/',sep='')


   sp.list=read.delim( "./QuasR_Aligned.txt",sep='\t')

   write.table(sp.list,'./sample_BAM.tmp',sep='\t',row.names=FALSE)

   cluObj=makeCluster(3)

   NOMEproj=qAlign(sampleFile='./sample_BAM.tmp', 
                genome="BSgenome.Celegans.UCSC.ce11", 
                paired="fr", 
                bisulfite="dir", 
                projectName="dSMF_N2", 
                clObj=cluObj)
   NOMEaln=as.data.frame(alignments(NOMEproj)[[1]])

samples=NOMEaln$SampleName



meth_gr=qMeth(NOMEproj, mode='allC',clObj=cluObj)

#dir.create(paste0(path,"/methylation_calls"))
saveRDS(meth_gr,paste0(path,'/methylation_calls/NomeAmpliconDM220216.rds'))

cO=3
PolII_in=call_context_methylation(meth_gr,cO,genome=Celegans)

hist(unlist(mcols(PolII_in[["CG"]])$V1),breaks=100)
hist(unlist(mcols(PolII_in[["GC"]])$V1),breaks=100)

saveRDS(PolII_in,paste0(path,"/methylation_calls/NomeAmpliconCEmeth.rds"))
meth_gr<-readRDS(paste0(path,"/methylation_calls/NomeAmpliconCEmeth.rds"))

#nci=rowSums(as.data.frame(elementMetadata(meth_gr))[,])>0
#meth_gr_f=meth_gr[nci]

designT=read.csv('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/finalChosenList_tss_dc_tissue.csv',
                stringsAsFactors=FALSE,header=TRUE)
design=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/amplicons_stranded.bed')
TSS=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/ampliconTSS.bed')
#most5p=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/most5p.bed')
#start(most5p)<-end(most5p)

#targeted
#designX=readRDS('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/XdcPrimers/tmp_design_test.rds')
#designA=readRDS('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/AdcPrimers/tmp_design_test.rds')
#design2=readRDS('/work2/gschub/arnaud/amplicon_sequencing/nomeseq/drosophila/tmp_design_PolII_V2.rds')
#design1= readRDS('/work2/gschub/arnaud/amplicon_sequencing/nomeseq/drosophila/DM_DE_val_ref_gr.rds')
#design3=readRDS('/work2/gschub/arnaud/nomeseq/WG/DM/analysis/primers/DM_DE_promoter_LPS_primers.rds')[[2]]
#design4=readRDS('/work2/gschub/arnaud/nomeseq/WG/DM/amplicon/tmp_design_Ecd.rds')
#reg=unlist(GRangesList(lapply(design4,'[[',2)))
# 	write.table(cbind(as.data.frame(reg)[,1:3]	,rep('',length(reg)),rep('',length(reg)),as.data.frame(reg)[,5]),'./rds/Ecd_primers_bed',sep='\t',quote=F,row.names=F,col.names=F)

#designA=unlist(GRangesList(lapply(designA,'[[',2)))

#design1=GRanges(seqnames(design1),ranges(design1),strand=strand(design1))
#designX=unlist(GRangesList( lapply(designX,'[[',2))) 
#amplicons=unique(c(designX,designA))
amplicons=design
amplicons=resize(amplicons,600,fix='center')
cluObj=makeCluster(3)

#call methylation
#TSS=import.bed('~/Documents/MeisterLab/dSMF/PromoterPrimerDesign/scripts/tss.bed')
#TSS=most5p

TSS2=resize(TSS,600,fix='center')
strand(TSS2)='+'

meth_gr <- qMeth(NOMEproj,mode="allC",TSS2)
saveRDS(meth_gr,paste0(path,'/Amplicon_raw_methylation_PolII_covered_amplicons.rds'))

#meth_gr <- qMeth(NOMEproj,mode="allC",TSS2)
#saveRDS(meth_gr,paste0(path,'/Amplicon_raw_methylation_PolII_covered_amplicons.rds'))


#meth_gr_f=readRDS('/work2/gschub/arnaud/nomeseq/WG/DM/amplicon/methylation_calls/NomeAmpliconDM250614.rds')

cO=10 #minimum coverage of site
PolII_in=call_context_methylation(meth_gr,cO,Celegans)

PolII_inm=unlist(GRangesList(PolII_in))
# saveRDS(PolII_inm,'/work2/gschub/arnaud/nomeseq/WG/DM/amplicon/methylation_calls/PolII_inhib_av10.rds')
saveRDS(PolII_inm,paste0(path,'PolII_inhib_av10.rds'))
PolII_inm=readRDS(paste0(path,'PolII_inhib_av10.rds'))


naix=apply(as.matrix(elementMetadata(PolII_inm[,1])),1,function(x){sum(is.na(x))==1})

#PolII_inm<-PolII_inm[!naix]
#saveRDS(PolII_inm,paste0(path,'PolII_inhib_av10.rds'))


Fl.s=c(grep('V1',colnames(elementMetadata(PolII_inm))),grep('V1',colnames(elementMetadata(PolII_inm))))
#Trp.s=c(grep('Trp',colnames(elementMetadata(PolII_inm))),grep('DMSO',colnames(elementMetadata(PolII_inm))))

pairs(1-as.matrix(elementMetadata(PolII_inm[!naix,Fl.s])),pch='.')
#pairs(1-as.matrix(elementMetadata(PolII_inm[!naix,Trp.s])),pch='.')
#cor(na.omit(1-as.matrix(elementMetadata(PolII_inm[!naix,Trp.s]))))
cor(na.omit(1-as.matrix(elementMetadata(PolII_inm[!naix,Fl.s]))))

library(gplots)
heatmap.2(cor(na.omit(1-as.matrix(elementMetadata(PolII_inm[!naix,Trp.s])))),trace='none',margins=c(10,10))
heatmap.2(cor(na.omit(1-as.matrix(elementMetadata(PolII_inm[!naix,Fl.s])))),trace='none',margins=c(10,10))
sp.names= colnames(elementMetadata(PolII_inm[!naix,1]))
outMet=sort(PolII_inm[!naix])

lapply(sl(sp.names),function(i){
   grToIgv(outMet,colnames(elementMetadata(outMet[,i])),paste(path,'/amplicon_DE_',sp.names[i],'.bed',sep=''))
})



#look hom many amplicons are covered

#amplicons= readRDS('/work2/gschub/arnaud/amplicon_sequencing/nomeseq/drosophila/DM_DE_val_ref_gr.rds')
#amplicons=resize(amplicons,600,fix='center')

table(countOverlaps(amplicons,PolII_inm[!naix,],ignore.strand=T)>0)

#FALSE  TRUE 
#11    85 

AmpOL<-findOverlaps(amplicons,PolII_inm[!naix,],ignore.strand=T)


pdf("AmpliconViews_AmpliconBiSeq.pdf",paper="a4",height=11,width=8)
for (i in 1:length(amplicons)) {
   b<-AmpliconViews(NOMEproj,amplicons[i],sampleNames="N2_DE_ampV001")
   try(plotAmpliconView(b),silent=TRUE)
}
dev.off()


a<-getCMethMatrix(NOMEproj,amplicons[1],"N2_DE_ampV001")

plotAmpliconMeth<-function(methMatrix,amplicon,TSS)  {
   numReads=dim(methMatrix)[1]
   plot(x=c(start(amplicon),end(amplicon)), y=c(1,numReads+1),type="n",
        xlab="position",ylab="read no.")
   methPos<-as.numeric(colnames(methMatrix))
   for (l in 1:numReads) {
      naidx<-is.na(methMatrix[l,])
      lineCol=ifelse(methMatrix[l,]==1,"black","gray")
      lineCol[naidx]<-NA
      for (p in 1:methPos) {
         lines(c(methPos[p],methPos[p]),c(l,l+1),col=lineCol[p])
      }
   }
   abline(TSS,col="red")
}
