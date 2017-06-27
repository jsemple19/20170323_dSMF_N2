library(QuasR)
library("BSgenome.Celegans.UCSC.ce11")

# collect citations for packages used
packageBib<-toBibtex(c(citation("QuasR"),citation("BSgenome.Celegans.UCSC.ce11")))

setwd("~/Documents/MeisterLab/sequencingData/20170323_dSMF_N2/dSMFseq_scripts")

source('./callAllCs.r') #to call all Cs
source('./useful_functionsV1.r') #load the ranges

#library("BSgenome.Ecoli.NCBI.20080805")  # delete this line??

# make clusters, needed to run in parallel
cluObj=makeCluster(3)

#setup directory structure this to desired location of your alignments
if (!dir.exists("../aln")){
   dir.create("../aln")
}
if (!dir.exists("../tmp")) {  #for trimmomatic quality trimmed reads
   dir.create("../tmp")
}
if (!dir.exists("../rds")) {
  dir.create("../rds")
}
if (!dir.exists("../plots")) {
  dir.create("../plots")
}
if (!dir.exists("../bed")) {
  dir.create("../bed")
}
path='../'
my.alignmentsDir=paste(path,'aln/',sep='')
tmp=paste(path,'tmp/',sep='')

#load the experiments
seqExp=read.table(paste0(path,'rawData/sampleList.txt'),sep='\t',header=T)

#create the QuasR Aln file
samples=as.data.frame(cbind(FileName1=paste(path,"rawData/",seqExp$FileName1,sep=''),
                            FileName2=paste(path,"rawData/",seqExp$FileName2,sep=''),
                            SampleName=as.character(seqExp$SampleName)))


###########################
#trim the low quality bases
###########################

for(i in sl(samples[,1])){
    #i=1
   spID=as.character(samples$SampleName[i])
   #clip the low quality bases #remove adapters
   system(paste(
      'java -jar /home/jenny/Trimmomatic-0.32/trimmomatic-0.32.jar PE ',
      samples$FileName1[i],' ', samples$FileName2[i], ' ',
      '../tmp/',samples$SampleName[i],'_forward_paired.fq.gz ',
      '../tmp/',samples$SampleName[i],'_forward_unpaired.fq.gz ',
      '../tmp/',samples$SampleName[i],'_reverse_paired.fq.gz ',
      '../tmp/',samples$SampleName[i],'_reverse_unpaired.fq.gz ',
      'ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36',
      sep='')
   )
}


AlnInput=as.data.frame(cbind(
   FileName1=paste(tmp,samples$SampleName,'_forward_paired.fq.gz',sep=''),
   FileName2=paste(tmp,samples$SampleName,'_reverse_paired.fq.gz',sep=''),
   SampleName=as.character(samples$SampleName)))

write.table(AlnInput,'QuasR_input.txt',quote=F,row.names=F,sep='\t')



###########################
#Align the full length fragments
###########################
QuasRdef='-k 2 --best --strata'

NOMEproj=qAlign(sampleFile=paste(path,"dSMFseq_scripts/QuasR_input.txt",sep=''),
                genome="BSgenome.Celegans.UCSC.ce11",
                #auxiliaryFile="/work2/gschub/Juliane/scripts/Altuna/auxFile.forQuasR.txt",
                aligner="Rbowtie",
                paired="fr",
                bisulfite="undir",
                projectName="dSMF_N2",
                alignmentsDir=my.alignmentsDir,
                clObj=cluObj ,
                alignmentParameter=paste('-e 150 -X 600 ',QuasRdef,sep=''),
                cacheDir = tmp)


#QC of the alignment

qQCReport(NOMEproj,'QC_QualityTrimmed_130616.pdf',clObj=cluObj)


alignments1=as.data.frame(alignments(NOMEproj)[[1]])




write.table(alignments1,'QuasR_Aligned.txt',quote=F,col.names=F,row.names=F,sep='\t',append=T)


