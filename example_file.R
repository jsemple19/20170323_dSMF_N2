
source('/work2/gschub/arnaud/Rscripts/useful_functionsV1.r') #load the ranges
source('/work2/gschub/arnaud/Rscripts/callAllCs.r') #to call all Cs


#############
#Functions
#############
##############################
#vectorize reads for plotting
##############################	

		VectorizeReads=function(target,mergedMat){ #vectorize reads from matrices in order to plot in base space
		 
# 			 mergedMat=matList[[i]][naidx,][h$rowInd,]
			GCmet=as.vector(t(mergedMat))
			GCstartP=rep(as.numeric(colnames(mergedMat)),nrow(mergedMat))
			GCreadID=unlist(lapply(seq(nrow(mergedMat)),function(i){rep(i,ncol(mergedMat))}))
			GCgr=GRanges(seqnames(target),IRanges(as.numeric(colnames(mergedMat)),as.numeric(colnames(mergedMat))))
			list(GCstartP,GCreadID,GCmet)
			}


		WGTiles<-function(genome=Dmelanogaster,spName){
			regs <- tileGenome(seqlengths(genome),tilewidth=1000000,cut.last.tile.in.chrom=TRUE)
			regsD <- cbind(as.data.frame(regs)[,1:3],sample=spName)
			regsD
# 	  	regsL <- split(regsD,1:nrow(regsD))
		}
	
	sorTSSbyReg <- function(regDF,sampleSheet,target_range,colBins=list(INRs=c(-16,4),DPEs=c(28,47),BREs=c(-36,-20),UPs=c(-58,-43))){
		#regDF:
		#sampleSheet: QuasR alignment file
			#the samples have to mention the treatment type by _NO_ _SS_ _DE_
		#target_range:
		#inMs|upMs|doMs: the collecting bins (relative to center)
	
			#######
			#collect the methylation data
			#######
			  projAll <- qAlign(sampleSheet,"BSgenome.Dmelanogaster.UCSC.dm3",paired="fr",bisulfite="dir")
			  proj <- projAll[alignments(projAll)$genome$SampleName==regDF[1,4]]
# 				regDF=regsL[[2]]
		
			  reg <- GRanges(seqnames = Rle(regDF[1,1]),ranges = IRanges(regDF[1,2], end = regDF[1,3]),seqlengths=seqlengths(Dmelanogaster))
			  regExp <- reg # expand the region at the end to catch all the pairs within the window
			  end(regExp) <- pmin(end(regExp)+2000,seqlengths(regExp)[as.character(seqnames(regExp))])

			  dm3_GCs_chr_XSV <- matchPattern("GC",Dmelanogaster[[as.character(seqnames(reg))]])
			  dm3_GCs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_GCs_chr_XSV), end = end(dm3_GCs_chr_XSV)),strand="+")
			  start(dm3_GCs_chr) <- start(dm3_GCs_chr)+1
			  dm3_GCs_chr_coord <- start(dm3_GCs_chr)

			  dm3_CGs_chr_XSV <- matchPattern("CG",Dmelanogaster[[as.character(seqnames(reg))]])
			  dm3_CGs_chr <- GRanges(seqnames = Rle(as.character(seqnames(reg))),ranges = IRanges(start(dm3_CGs_chr_XSV), end = end(dm3_CGs_chr_XSV)),strand="+")
			  end(dm3_CGs_chr) <- end(dm3_CGs_chr)-1
			  dm3_CGs_chr_coord <- start(dm3_CGs_chr)

			  methAln <- qMeth(proj,regExp,mode="allC",reportLevel="alignment")[[1]]
				
				if(length(  methAln[[1]])!=0){ #checks that the bin hase some reads (important for amplicon)
					
					  if(length(grep("_NO_",regDF[1,4])) > 0){
						#GC
						methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]+1
						selCs <- methAln$Cid %in% dm3_GCs_chr_coord
					  }else if(length(grep("_SS_",regDF[1,4])) > 0){
						#CG
						methAln$Cid[methAln$strand=="-"] <- methAln$Cid[methAln$strand=="-"]-1
						selCs <- methAln$Cid %in% dm3_CGs_chr_coord
					  }else if(length(grep("_DE_",regDF[1,4])) > 0){
						#GC&CG simultaneously. strand collapse is more complicated. prioritize GC in GCG minus strand collapse
						plusStrandSel <- methAln$strand=="+"
						SelOvGC_plus <- ((methAln$Cid %in% dm3_GCs_chr_coord) & plusStrandSel)
						SelOvGC_minus <- ((methAln$Cid+1) %in% dm3_GCs_chr_coord) & !plusStrandSel
						SelOvCG_plus <- ((methAln$Cid %in% dm3_CGs_chr_coord) & plusStrandSel)
						SelOvCG_minus <- ((methAln$Cid-1) %in% dm3_CGs_chr_coord) & !plusStrandSel
						methAln$Cid[SelOvGC_minus] <- methAln$Cid[SelOvGC_minus]+1
						methAln$Cid[!SelOvGC_minus & SelOvCG_minus] <- methAln$Cid[!SelOvGC_minus & SelOvCG_minus]-1
						selCs <- SelOvGC_plus | SelOvGC_minus | SelOvCG_plus | SelOvCG_minus
					  }else{stop("Sample name error. No _McvPI_ , _SsssI_ or _McvPISsssI_")}
	
	
					#######
					#Sort the reads
					######
			
					regT=GRanges( regDF[1,1],IRanges(regDF[1,2],regDF[1,3]))
					TFs=subsetByOverlaps(target_range,regT,type='within')
			
					sRs=	sortTSSreads_vect(methAln,selCs,TFs,colBins)
					sRs
						}else{c()}
				}
 
 
 
 
 
 library("BSgenome.Dmelanogaster.UCSC.dm3")
# source('/work2/gschub/arnaud/Rscripts/Single_molecule_manipulation _function.r')
source('/work2/gschub/arnaud/Rscripts/useful_functionsV1.r')
# library(BSgenome)

 
 
	sortTSSreads_vect<-function(methAln,selCs,st,colBins){
				#Daughter function called by 	sorTFbyReg 
				#methAln: qMeth object
				#selCs: Cs in the correct context (to be used for sorting)
				#st: the motifs cooridnate GR
				#inMs|upMs|doMs: the collecting bins (relative to center)
			
		 #create  Crange	
			Crange=IRanges(methAln[[2]][selCs],methAln[[2]][selCs])
			#create collecting intervals
			#st=refGR[hspids][5]
			st=resize(st,1,fix='center')
			midP=start(st)
		
# 			INRMP=GRanges(seqnames(st),IRanges(ifelse(as.logical(strand(st) =='+'),midP+INRs[1],midP-INRs[2]),ifelse(strand(st)=='+',midP+INRs[2],midP-INRs[1])))
			#create the ranges of all bins	
			binMP=lapply(sl(colBins),function(b){
				GRanges(seqnames(st),IRanges(ifelse(as.logical(strand(st) =='+'),midP+colBins[[b]][1],midP-colBins[[b]][2]),ifelse(strand(st)=='+',midP+colBins[[b]][2],midP-colBins[[b]][1])))
			})
			#find the C overlaps within bins			
			binO=lapply(sl(colBins),function(b){
				as.matrix(findOverlaps(ranges(binMP[[b]]),Crange)  )
			})
			
			
		   #compute the methylation vectors
		   binM=lapply(sl(colBins),function(b){
				methAln[[4]][selCs][binO[[b]][,2]]
			})
		
			#compute the IDs
			binID=lapply(sl(colBins),function(b){
				paste(names(st)[binO[[b]][,1]],methAln[[1]][selCs][binO[[b]][,2]],sep='_')
			})
		
			
			#group Cs per region/per_read
			g.binM=lapply(sl(colBins),function(b){
				round(tapply(   binM[[b]],binID[[b]],mean))
			})
			u.binID=lapply(sl(colBins),function(b){
				sort(unique(binID[[b]]))
			})
		
			
			#intersect IDS
			Uids=unique(unlist(u.binID))
			o.mat=do.call(cbind,lapply(sl(colBins),function(b){
					Uids%in%u.binID[[b]]
					}))
			sh.ids=Uids[rowSums(o.mat)==ncol(o.mat)]
			
			if(length(sh.ids)>0){ #check if there are reads spanning the whole region
				#binary met vectors
				sh.ninM=lapply(sl(colBins),function(b){
					g.binM[[b]][sh.ids]
				})
				met.mat=data.frame(do.call(cbind,sh.ninM))
		
	# 			pattern=paste(as.character(s.inM),as.character(s.outM),sep='')
				pattern=apply(met.mat,1,function(x){paste(as.character(x),sep='',collapse='')})
			
				st.id=string.split(sh.ids,'_',1)
				read.id=string.split(sh.ids,'_',2)
	
				#SORTED READ LISTS
				sR=split(read.id,paste(st.id,pattern,sep='_'))
				sRl=split(sR,string.split(names(sR),'_',1))
	
				sRl
			}else{c()}
		}
  			



#####################
#############
#example call
#############
################




		sampleSheet='/work2/gschub/arnaud/HTS/bis-seq/NOME-seq/DM/meta_aln/QuasR_AlignedSubset.txt'
		dm3_refseq_TSS.cs=readRDS('/work2/gschub/arnaud/nomeseq/WG/DM/analysis/tmp_refseq_TSS_S2_CAGE.rds')

		
		#Sort reads

		for(i in (1:length(spNames))[10:13]){
			#i=10
			regsD=WGTiles(Dmelanogaster,spNames[i])
			regsD=regsD[regsD[,1]%in%c('chr4','chrX','chr2L','chr2R','chr3L','chr3R'),]
			regsL <- split(regsD,1:nrow(regsD))
		
			#set the reference object
			refGR=resize(dm3_refseq_TSS.cs,1,fix='center')
		
			sTSS.l=mclapply(sl(regsL),function(reg.i){
				regDF=regsL[[reg.i]]
				sTSS=sorTSSbyReg(regDF,sampleSheet,	refGR,colBin)
				},mc.cores=15)
			sTSS.lu=unlist(sTSS.l,recursive=F)
			 saveRDS(list(refGR,sTSS.lu),paste('./rds/tmp_sR_TSS_',spNames[[i]],'_WG.rds',sep=''))
			}
		
		#esclude simple treatment
		spNames=spNames[grep('_DE_',spNames)]
		
		sR_TSS=mclapply(sl(spNames),function(i){
		
			readRDS(paste('./rds/tmp_sR_TSS_',spNames[[i]],'_WG.rds',sep=''))[[2]]
			},mc.cores=5)
		names(sR_TSS)=spNames
		

		refGR=readRDS(paste('./rds/tmp_sR_TSS_',spNames[[2]],'_WG.rds',sep=''))[[1]]
 		sortedReadsList=sR_TSS


			#Sort reads
# 			#create a list of count matrices
			countMats=lapply(sl(spNames),function(i){
				countMat=matrix(NA,nrow=length(refGR),2^4)
				allPos=expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
				patternStrings=names(table(apply(allPos,1,function(x){(paste(as.character((x)),collapse=''))})))
				colnames(countMat)=patternStrings
				countMat
				})
			names(countMats)=spNames
	
			#populate the matrices
			for (spi in sl(spNames)){
			print(spi)
				readIdList=sortedReadsList[[spi]]
				
				for (i in as.character(sl(refGR))){
 					if(i %in%  names(readIdList)){
							IDs=readIdList[[i]]
							counts=unlist(lapply(IDs,length))
							countMats[[spi]][as.numeric(i),]=0
							countMats[[spi]][as.numeric(i),string.split(names(IDs),'_',2)]=counts
 							
 							}}
 							}
				
 	
 	
 
 	
 	saveRDS(countMats,'./rds/tmp_WG_refseq_countMats.rds')
# 	countMats=readRDS('./rds/tmp_WG_TSS_countMats.rds')
 	countMats=readRDS('./rds/tmp_WG_refseq_countMats.rds')


		#get a frequency matrix (individual states)
		cutoff=20
		freqMat=mclapply(sl(spNames),function(spi){
				totC=1+apply(countMats[[spi]],1,function(x){sum(x,na.rm=T)})
	
				apply((countMats[[spi]]),2,function(x){
					fid=totC>cutoff
					y=(x/totC)*100
					y[!fid]<-NA
					y
				})},mc.cores=10)
		
`	names(freqMat)=spNames
