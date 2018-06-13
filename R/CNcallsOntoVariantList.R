#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
variantFile=args[1]
ASCATfile=args[2]
ploidyFile=args[3]

library(GenomicRanges)

# add copy number calls from ASCAT to variant sheet
addCN = function(seg,segGrange,variantGrange,ploidy,thresh=1.3)
	{
	# Granges
	#segGrange = as(paste0("chr",seg[,2],":",seg[,3],"-",seg[,4]),"GRanges")
	#variantGrange = as(paste0("chr",chrom,":",start,"-",end),"GRanges")
	# overlap
	overlap = findOverlaps(segGrange,variantGrange)
	seg = seg[overlap@from,,drop=FALSE]
	if(nrow(seg)==0) 
		{
		out = rep(NA,5)
		names(out) = c("tumCNtot","tumCNb","ploidy","log2ploidyTum","CNcalls")
		return(out)
		}
	tag = NULL
	# sep CNs
	CNab = cbind(seg[,7]-seg[,8],seg[,8])
	# LOH
	LOHcheck = apply(CNab,MARGIN=1,FUN=function(x) sum(x==0)==1)
	if(LOHcheck) tag = "LOH"
	# Hom del
	if(any(seg[,7]==0)) tag = "Hom del"
	# amplified
	if(all(log2(seg[,7]/ploidy)>thresh)) tag = paste0(c(tag,"CN gain"),collapse=";")
	# deleted
	if(any(log2(seg[,7]/ploidy)<c(-thresh))) tag = paste0(c(tag,"CN loss"),collapse=";")
	if(length(tag)==0) tag=NA
	# out object
	out = c(paste(seg[,7],collapse=";"),
		paste(seg[,8],collapse=";"),
		ploidy,
		paste(log2(seg[,7]/ploidy),collapse=";"),
		tag)
	names(out) = c("tumCNtot","tumCNb","ploidy","log2ploidyTum","CNcalls")
	return(out)
	}

# ASCAT CN
CN = read.table(ASCATfile,as.is=TRUE,sep="\t")
CNgRange = as(paste0("chr",CN[,2],":",CN[,3],"-",CN[,4]),"GRanges")
# variants sheet
variants = read.table(variantFile,quote="",sep="\t",comment.char="@",as.is=TRUE,head=TRUE)
variantGrange = as(paste0("chr",variants[,"Chrom"],":",variants[,"Start_Position"],"-",variants[,"End_Position"]),"GRanges")
# ploidy estimates
ploidy = read.csv(ploidyFile)
# run CN
toRun = 1:nrow(variants)
CNcalls = sapply(toRun, FUN=function(x) {print(x);
	addCN(seg=CN[which(CN[,1]==variants[x,"Sample"]),,drop=FALSE],
		segGrange = CNgRange[which(CN[,1]==variants[x,"Sample"])],
		variantGrange = variantGrange[x],
		ploidy=ploidy$Ploidy[which(ploidy$Sample==variants[x,"Sample"])]
		)})

variants = cbind(variants,t(CNcalls))
# get out file name
fileName = rev(strsplit(variantFile,split="/")[[1]])[1]
fileName = strsplit(fileName,split="[.]")[[1]][1]
fileNameOut = paste0(fileName,"-withCNcalls.txt")
outFile = gsub(fileName,fileNameOut,variantFile)
# write out file
write.table(variants,file="/media/grp11/CI_Pathology_Steele/projects/1_undiff/results/WGS/annot/combined/v2-16-03-2017/combined_sanger_vep_annot_no_first_variant_with_CADD-filteredNonCoding-withPriority-withCNcalls.txt",sep="\t",quote=FALSE,row.names=FALSE,na="")
