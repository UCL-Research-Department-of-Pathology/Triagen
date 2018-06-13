#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test arguments: if not, return an error
if (!length(args)%in%c(4,16)) 
	{
	stop("Incorrect number of arguments", call.=FALSE)
	} else {
		if (length(args)==1) 
			{
			scriptDir=dirname(sys.frame(1)$ofile)
			setwd(scriptDir)
			args[2:4] = c("../data/knownMuts/genie_known_muts.bed",
				"../data/knownMuts/civic_variants_03022017.tsv",
				"../data/knownMuts/Sanger_drivers.csv")
			} else if(length(args!=16)) {
			# default column headings
			args[5:18] = c(
				"Chrom",				# chromosome
				"Start_Position",			# start
				"End_Position",				# end
				"No.of.Occurances.in.TCGA.ICGC",	# dbCountCol 
				"ExAC_AF",				# exacCountCol 
				"Hugo_Symbol",				# geneCol 
				"HGVSc",				# variantCol 
				"IMPACT",				# impactCol 
				"CADD_PHRED",				# caddCol
				"CLIN_SIG",				# clinsigCol
				"Variant_Classification",		# variantClass
				"Sample",				# patientCol
				"Ref",					# reference
				"Alt"					# alternative
				)
			}
		}
	}


EXOME <<- 0
EXAC <<- 0
GENIE <<- 0

TCGAZERO <<- 0
TCGAZEROOBSULC <<- 0
TCGAZEROOBSLC <<- 0
TCGAZEROPREDUNRELIABLE <<- 0
TCGAZEROPREDULC <<- 0
TCGAZEROPREDLC <<- 0

TCGALOW <<- 0
TCGALOWOBSULC <<- 0
TCGALOWOBSLC <<- 0
TCGALOWPREDUNRELIABLE <<- 0
TCGALOWPREDULC <<- 0
TCGALOWPREDLC <<- 0

TCGAHIGH <<- 0
PATHUNRELIABLE <<- 0
INDEX <<- 0

classifyVariant = function(chrom, 	# chromosome
			posStart, 	# start position
			posEnd, 	# end position
			ref,		# ref seq
			alt,		# alt seq
			cDNA,		# cDNA for variant
			nPatients, 	# count of variant in patients
			tcgaCount=NULL, # count of variant in TCGA and ICGC
			exacCount, 	# count of variant in exac
			impact, 	# VEP impact
			caddScore, 	# CADD score
			genie,		# genie hotspot mutations
			clinsig,	# clinical significance
			variantClass,	# variant classification
			civic,		# civic databse
			sanger,		# sanger database
			unidirectionalFilter, # unidirectionalFilter
			germlineFilter, germlineThresh=0.015)	# germlineFilter
	{
	# unidirectional
	if(unidirectionalFilter==0)
		{
		return(c("Unreliable","Unidirectional"))
		}
	# germline filter
	if(germlineFilter>germlineThresh)
		{
		return(c("Unreliable","Germline"))
		}
	# silent mutations for driver analysis
	if(variantClass=="Silent") return(c("silent","silent"))
	# exac
	if(!is.na(exacCount)) 
		{
		EXAC <<- EXAC+1
		return(c("Unreliable","ExAC"))
		}
	# genie
	variant = data.frame(chrom=chrom,start=as.numeric(posStart),end=as.numeric(posEnd))
	variant = as(variant,"GRanges")
	overlaps = findOverlaps(genie,variant)
	if(length(overlaps)>0) 
		{
		GENIE <<- GENIE+1
		return(c("High_confidence","Genie"))
		}
	# civic
	matchIndex = which(paste0(civic$chromosome)==paste0(chrom)&
			civic$start==posStart&
			civic$stop==posEnd&
			paste0(civic$reference_bases)==paste0(ref)&
			paste0(civic$variant_bases)==paste0(alt))
	if(length(matchIndex)>0) return(c("High_confidence","CIVIC"))
	# sanger
	matchIndex = which(paste0(sanger$Chr)==paste0(chrom)&
			sanger$Posn==posStart&
			paste0(sanger$cDNA)==paste0(cDNA))
	if(length(matchIndex)>0) return(c("High_confidence","Sanger"))
	# TCGA/ICGC
	if(!is.null(tcgaCount))
		{
		if(tcgaCount==0)
			{
			TCGAZERO <<- TCGAZERO + 1 
			category = triageCurrentObs(nSamples,nPatients)
			if(category=="ultra_low_confidence")
				{
				TCGAZEROOBSULC <<- TCGAZEROOBSULC + 1 
				} else {
				TCGAZEROOBSLC <<- TCGAZEROOBSLC + 1
				}
			category = triagePathPred(caddScore,impact,clinsig,category)
			if(category=="Unreliable") 
				{
				TCGAZEROPREDUNRELIABLE <<- TCGAZEROPREDUNRELIABLE + 1 
				return(category)
				} else if(category=="ultra_low_confidence") {
				TCGAZEROPREDULC <<- TCGAZEROPREDULC + 1 
				} else {
				TCGAZEROPREDLC <<- TCGAZEROPREDLC + 1 
				}
			return(c(paste0("Novel_",category),"TCGA/CADD"))
			} else if(tcgaCount<3) {
			TCGALOW <<- TCGALOW + 1
			category = triageCurrentObs(nSamples,nPatients)
			if(category=="ultra_low_confidence")
				{
				TCGALOWOBSULC <<- TCGALOWOBSULC + 1 
				} else {
				TCGALOWOBSLC <<- TCGALOWOBSLC + 1
				}
			category = triagePathPred(caddScore,impact,clinsig,category)
			return(c(paste0("Observed_",category),"TCGA/CADD"))
			} else {
			TCGAHIGH <<- TCGAHIGH + 1
			return(c("Medium_confidence","TCGA"))
			}
		} else {
		category = triageCurrentObs(nSamples,nPatients)
		category = triagePathPred(caddScore,impact,clinsig,category)
		return(c(category,"CADD"))
		}
	# something has gone wrong if reach here
	print("error")
	return("error")
	}

triageCurrentObs = function(nSamples,nPatients)
	{
	# more than one sample, more than one patient
	if(nPatients>1) return("low_confidence")
	# more than one sample, all in same patient
	return("ultra_low_confidence")
	}

triagePathPred = function(CADD,impact,clinsig,designation)
	{
	lowFlag = grepl("ultra",designation)	
	# increase priority of high impact 
	if(lowFlag&((CADD>30&!is.na(CADD))|impact=="HIGH"|clinsig=="pathogenic")) return(gsub("ultra_","",designation))
	# decrease priority of low impact 
	if(!lowFlag&((CADD<20&!is.na(CADD))&impact%in%c("LOW","MODIFIER"))) return(paste0("ultra_",designation))
	return(designation)
	}

# column names for columns of interest
chromCol = args[5]
posStartCol = args[6]
posEndCol = args[7]
tcgaCountCol = args[8]
exacCountCol = args[9]
geneCol = args[10]
variantCol = args[11] 
impactCol = args[12]
caddCol = args[13]
clinsigCol = args[14]
variantClassCol = args[15]
patientCol = args[16]
refCol = args[17]
altCol = args[18]
# data files
dataFile = args[1]
genieFile = args[2]
civicFile = args[3]
sangerFile = args[4]
# load data
print("read data")
data = read.table(dataFile,head=TRUE,sep="\t",comment.char="@",quote="",as.is=TRUE)
# load genie
print("read genie")
genie = read.csv(genieFile,head=FALSE,skip=10)
# convert genie to 1-based
genie[,2] = genie[,2]+1
colnames(genie)[1:3] = c("chrom","start","end")
library(GenomicFeatures)
genie = as(genie,"GRanges")
# load civic
print("read civic")
civic = read.table(civicFile,head=TRUE,sep="\t",comment.char="@",quote="")
civic = civic[which(civic$reference_bases!=""),]
# load sanger
print("read sanger")
sanger = read.csv(sangerFile,head=TRUE)
# get unidirectional filter
print("set unidirectional filter")
getMinStrand = function(data)
	{
	if(data$Type=="Sub")
		{
		# substitutions - minimum of alternative
		alt = data[,"Alt"]
		return(min(data[,paste0("F",alt,"Z.Tum")],data[,paste0("R",alt,"Z.Tum")]))
		} else {
		# indels - minimum of unique calls (Pindel and BWA)
		return(min(c(data$PU.Tum,data$NU.Tum)))
		}
	}
unidirectionalFlag = sapply(1:nrow(data),FUN=function(x) getMinStrand(data[x,,drop=FALSE]))
data = cbind(data,unidirectionalFlag)
uniCol = "unidirectionalFlag"
# get germline filter
print("set germline filter")
getGermline = function(data,infoMut=c("PU.Norm","NU.Norm"),infoAll=c("PR.Norm","NR.Norm"))
	{
	if(data$Type=="Sub")
		{
		# substitutions - normal alt / normal alt & ref
		alt = data[,"Alt"]
	        ref = data[,"Ref"]
		altN = sum(as.numeric(data[,paste0(c("F","R"),alt,"Z.Norm")]))
		refN = sum(as.numeric(data[,paste0(c("F","R"),ref,"Z.Norm")]))
		return(altN/(altN+refN))
		} else {
		# indels - normal mutant / normal total
		return(sum(as.numeric(unlist(data[,infoMut])))/sum(as.numeric(unlist(data[,infoAll]))))
		}
	}
germlineFlag = sapply(1:nrow(data),FUN=function(x) getGermline(data[x,,drop=FALSE]))
data = cbind(data,germlineFlag)
germCol = "germlineFlag"







# function to count the number of patients with the variant before classifying
prioritiseVariant = function(chrom,posStart,posEnd,
			gene,variant,
			tcgaCount=NULL,
			exacCount,
			impact,
			cadd,
			clinsig,
			variantClass,
			ref,
			alt,
			unidirectional,
			germline)
	{
INDEX <<- INDEX+1
	print(paste0(INDEX,". ",gene,":",variant,",",chrom,":",posStart,"-",posEnd))
	# subset to this gene and variant
	subData = data[which(data[,geneCol]==gene&data[,variantCol]==variant),]
	# count number of patients with this variant
	patientCount = length(unique(subData[,patientCol]))
	# classify variant priority
	priority = classifyVariant(
		  chrom=chrom,
		  posStart=posStart,
		  posEnd=posEnd,
		  cDNA=variant,
		  nPatients=patientCount,
                  tcgaCount=tcgaCount,
                  exacCount=exacCount,
		  impact=impact,
		  caddScore=cadd,
		  genie=genie,
		  clinsig=clinsig,
		  variantClass=variantClass,
		  sanger=sanger,
		  civic=civic,
		  ref=ref,
		  alt=alt,
		  unidirectionalFilter=unidirectional,
		  germlineFilter=germline)
	return(c(priority,patientCount))
	}

toRun = 1:nrow(data)
# perform classification
Priority = mapply(FUN=prioritiseVariant,
		    chrom=data[toRun,chromCol],
		    posStart=data[toRun,posStartCol],
		    posEnd=data[toRun,posEndCol],
                    gene=data[toRun,geneCol],
                    variant=data[toRun,variantCol],
                    tcgaCount=data[toRun,tcgaCountCol],
                    exacCount=data[toRun,exacCountCol],
		    impact=data[toRun,impactCol],
		    cadd=data[toRun,caddCol],
		    clinsig=data[toRun,clinsigCol],
		    variantClass=data[toRun,variantClassCol],
		    ref=data[toRun,refCol],
		    alt=data[toRun,altCol],
		    unidirectional=data[toRun,uniCol],
		    germline=data[toRun,germCol]
                    )
rownames(Priority) = c("priority","reason","patientCount")
priorityTable = table(Priority["priority",])
reasonTable = table(Priority["reason",])

# combine data and results
data = cbind(data,t(Priority))

# get out file name
fileEnding = rev(strsplit(dataFile,split="[.]")[[1]])[1]
outFileSave = rev(strsplit(dataFile,"/")[[1]])[1]
outFile = gsub(outFileSave,"",dataFile)
# increment version
chars = strsplit(outFileSave,split="")[[1]]
if(chars[1]=="v") 
	{
	chars[2] = as.numeric(chars[2])+1
	outFileEnd = paste0(chars,collapse="") 
	} else {
	outFileEnd = outFileSave
	}
outFileEnd = strsplit(outFileEnd,"[.]")[[1]][1]
# write results
write.table(data,
	file=paste0(outFile,outFileEnd,"-withPriority.txt"),
	quote=FALSE,row.names=FALSE,sep="\t")
write.table(data[which(Priority["priority",]!="Unreliable"),],
	file=paste0(outFile,outFileEnd,"-withPriority-noUnreliable.txt"),
	quote=FALSE,row.names=FALSE,sep="\t")













