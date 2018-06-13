#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
expectedNargs = 29
nArgs = length(args)
# test arguments: if not, return an error
if (!nArgs%in%c(1,6,expectedNargs)) 
	{
	stop("Incorrect number of arguments", call.=FALSE)
	} else {
		if(nArgs==1) 
			{
			#scriptDir=dirname(sys.frame(1)$ofile)
#scriptDir <- getSrcDirectory(function(dummy) {dummy})
#print(scriptDir)
			#setwd(scriptDir)
			args[2:5] = c("~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/genie_mutations_2.0.0_2017-11-27.txt",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/civic_variants_2017-11-01.txt",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/Sanger_drivers.csv",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/mskcc_hotspots_2017-11-27.txt",
				"~/Dropbox/PostDoc/Rlibs/Triagen/data/knownMuts/encode_blacklist_2017-11-27.bed.gz")
			} 
		if(nArgs!=expectedNargs) 
			{
			# default column headings
			args[7:expectedNargs] = c(
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
				"Alt",					# alternative
				"Type",					# variant type
				"Sub",					# substition name
				TRUE,					# bool: run unidirectional filter
				FALSE,					# bool: run germline filter
				TRUE,					# bool: run CADD filter
				"sanger",				# germline method
				"n_ref_count",				# normal reference count
				"n_alt_count",				# normal alternative count
				TRUE					# doParallel
				)
			}
		}



# function to prioritise variants based on a number of factors
classifyVariant = function(
			nPatients, 	# count of variant in patients
			tcgaCount=NULL, # count of variant in TCGA and ICGC
			exacCount, 	# count of variant in exac
			impact, 	# VEP impact
			caddScore, 	# CADD score
			clinsig,	# clinical significance
			variantClass,	# variant classification
			genie,		# genie hotspot mutations
			civic,		# civic databse
			sanger,		# sanger database
			mskcc,		# memorial sloan kettering database
			blacklist,	# encode blacklist regions 
			unidirectionalFilter, # unidirectionalFilter
			germlineFilter, germlineThresh=0.015, # germlineFilter
			doUni=TRUE,doGermline=FALSE,doCADD=TRUE, # running options
			exacThresh=0.0004) # exac threshold
	{
	# civic
	if(civic) return(c("High_confidence","CIVIC"))
	# sanger
	if(sanger) return(c("High_confidence","Sanger"))
	# genie
	if(genie) return(c("High_confidence","Genie"))
	# mskcc
	if(mskcc) return(c("High_confidence","MSKCC"))
	# blacklist
	if(blacklist) return(c("Unreliable","Encode_blacklist"))
	# unidirectional
	if(doUni)
		{
		if(unidirectionalFilter==0)
			{
			return(c("Unreliable","Unidirectional"))
			}
		}
	# germline filter
	if(doGermline)
		{
		if(!is.na(germlineFilter))
			{
			if(germlineFilter>germlineThresh)
				{
				return(c("Unreliable","Germline"))
				}
			}
		}
	# silent mutations for driver analysis
	if(variantClass=="Silent") return(c("silent","silent"))
	# exac
	if(!is.na(exacCount)&exacCount>exacThresh) 
		{
		return(c("Unreliable","ExAC"))
		}
	# TCGA/ICGC
	if(!is.null(tcgaCount))
		{
		if(tcgaCount==0)
			{
			category = triageCurrentObs(nSamples,nPatients)
			if(doCADD) category = triagePathPred(caddScore,impact,clinsig,category)
			return(c(paste0("Novel_",category),"TCGA/CADD"))
			} else {
			category = triageCurrentObs(nSamples,nPatients)
			if(doCADD) category = triagePathPred(caddScore,impact,clinsig,category)
			return(c(paste0("Observed_",category),"TCGA/CADD"))
			}
		} else {
		category = triageCurrentObs(nSamples,nPatients)
		if(doCADD) category = triagePathPred(caddScore,impact,clinsig,category)
		return(c(category,"CADD"))
		}
	# something has gone wrong if reach here
	print("error")
	return("error")
	}

# function to triage base on numbers of observations in samples/patients
triageCurrentObs = function(nSamples,nPatients)
	{
	# more than one sample, more than one patient
	if(nPatients>1) return("low_confidence")
	# more than one sample, all in same patient
	return("ultra_low_confidence")
	}

# function to triage based on predicted pathogenicity
triagePathPred = function(CADD,impact,clinsig,designation)
	{
	lowFlag = grepl("ultra",designation)	
	# increase priority of high impact 
	if(lowFlag&((CADD>20&!is.na(CADD))|impact%in%c("HIGH","MODERATE")|(!is.na(clinsig)&clinsig=="pathogenic")))
		{
		return(gsub("ultra_","",designation))
		}
	# decrease priority of low impact 
	if(!lowFlag&((CADD<20&!is.na(CADD))&impact%in%c("LOW","MODIFIER")))
		{
		return(paste0("ultra_",designation))
		}
	return(designation)
	}


# function to check validity of columns
checkCol = function(data,column)
	{
	check = column%in%colnames(data)
	if(!check) stop(paste0("Column ",column," missing from data"))
	if(all(is.na(data[,column]))) print(paste0("Warning: all data values are NA in column ",column))
	}

# function to check files
checkFile = function(fileName)
	{
	# check if file exists
	if(!file.exists(fileName)) stop(paste0("File ",fileName," does not exist"))
	# check if file is empty
	info = file.info(fileName)
	if(info$size==0) stop(paste0("File ",fileName," is empty"))
	}

# get out filename
getOutFile = function(filename)
	{
	# get out file name
	print("set file name")
	fileEnding = rev(strsplit(filename,split="[.]")[[1]])[1]
	outFileSave = rev(strsplit(filename,"/")[[1]])[1]
	outFile = gsub(outFileSave,"",filename)
	# increment version
	chars = strsplit(outFileSave,split="")[[1]]
	if(chars[1]=="v") 
		{
		chars[2] = as.numeric(chars[2])+1
		outFileEnd = paste0(chars,collapse="") 
		} else {
		outFileEnd = outFileSave
		}
	# return filename
	return(paste0(outFile,outFileEnd,"-withPriority.txt"))
	}

# function to get germline filter
getGermlineSanger = function(data,infoMut=c("PU.Norm","NU.Norm"),infoAll=c("PR.Norm","NR.Norm"),
		variantTypeCol,altCol,refCol,subName)
	{
	if(data[,variantTypeCol]==subName)
		{
		# substitutions - normal alt / normal alt & ref
		alt = data[,altCol]
	        ref = data[,refCol]
		altN = sum(as.numeric(data[,paste0(c("F","R"),alt,"Z.Norm")]))
		refN = sum(as.numeric(data[,paste0(c("F","R"),ref,"Z.Norm")]))
		return(altN/(altN+refN))
		} else {
		# indels - normal mutant / normal total
		return(sum(as.numeric(unlist(data[,infoMut])))/sum(as.numeric(unlist(data[,infoAll]))))
		}
	}


# function to get unidirectional filter
getMinStrand = function(data,variantTypeCol,altCol,subName)
	{
	if(data[,variantTypeCol]==subName)
		{
		# substitutions - minimum of alternative
		alt = data[,altCol]
		return(min(data[,paste0("F",alt,"Z.Tum")],data[,paste0("R",alt,"Z.Tum")]))
		} else {
		# indels - minimum of unique calls (Pindel and BWA)
		return(min(c(data$PU.Tum,data$NU.Tum)))
		}
}

# function to set up prioritise
setupPrioritise = function(dataFile,genieFile,civicFile,sangerFile,mskccFile,blacklistFile,chromCol,
		posStartCol,posEndCol,tcgaCountCol,exacCountCol,geneCol,
		variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni=TRUE,doGermline=FALSE,doCADD=TRUE,
		germlineMeth="sanger",germlineRefCountCol=NA,germlineAltCountCol=NA)
	{	
	# check input files
	sapply(c(dataFile,
		genieFile,
		civicFile,
		sangerFile,
		mskccFile,
		blacklistFile),FUN=checkFile)
	# load data
	print("read data")
	data = read.table(dataFile,head=TRUE,sep="\t",comment.char="@",quote="",as.is=TRUE)
	# check data columns
	print("check data")
	sapply(c(chromCol,
		posStartCol,
		posEndCol,
		tcgaCountCol,
		exacCountCol,
		geneCol,
		variantCol,
		impactCol,
		clinsigCol,
		variantClassCol,
		refCol,
		altCol,
		variantTypeCol),FUN=function(x) checkCol(data,x))
	if(doCADD) checkCol(data,caddCol) 
	# load genie
	print("read genie")
	genie = read.csv(genieFile,head=FALSE,skip=10)
	print("genie check")
	dataCheck = paste(data[,chromCol],data[,posStartCol],data[,posEndCol],data[,refCol],data[,altCol],sep=":")
	genieCheck = paste(genie$Chromosome,genie$Start_Position,genie$End_Position,genie$Reference_Allele,genie$Tumor_Seq_Allele2,sep=":")
	genieMatch = dataCheck%in%genieCheck	
	# load civic
	print("read civic")
	civic = read.table(civicFile,head=TRUE,sep="\t",comment.char="@",quote="",as.is=TRUE)
	civic = civic[which(civic$reference_bases!=""&civic$variant_bases!=""),]
	civic$reference_bases[which(civic$reference_bases=="")] = "-"
	civic$variant_bases[which(civic$variant_bases=="")] = "-"
	print("civic check")
	civicCheck = paste(civic$chromosome,civic$start,civic$stop,civic$reference_bases,civic$variant_bases,sep=":")
	civicMatch = dataCheck%in%civicCheck
	# load mskcc
	print("read mskcc")
	mskcc = read.table(mskccFile,sep="\t",head=TRUE)
	print("mskcc check")
	mskccCheck = paste(mskcc$Chromosome,mskcc$Start_Position,mskcc$End_Position,mskcc$Reference_Allele,mskcc$Tumor_Seq_Allele2,sep=":")
  mskccMatch = dataCheck%in%mskccCheck
	# load sanger
	print("read sanger")
	sanger = read.csv(sangerFile,head=TRUE)
	print("sanger check")
	dataCheck = paste(data[,chromCol],data[,posStartCol],data[,variantCol],sep=":")
	sangerCheck = paste(sanger$Chr,sanger$Pos,sanger$cDNA,sep=":")
	sangerMatch = dataCheck%in%sangerCheck
	# load blacklist
	print("read blacklist")
	blacklist = read.table(blacklistFile,sep="\t",head=FALSE)
	colnames(blacklist)[1:3] = c("chrom","start","end")
	library(GenomicFeatures)
	blacklist = as(blacklist,"GRanges")
	print("blacklist check")
	variant = data.frame(chrom=paste0("chr",data[,chromCol]),start=as.numeric(data[,posStartCol]),end=as.numeric(data[,posEndCol]))
	variant = as(variant,"GRanges")
	overlaps = findOverlaps(blacklist,variant)
	blacklistMatch = sapply(1:nrow(data),FUN=function(x) x%in%overlaps@to)
	# get unidirectional filter
	print("set unidirectional filter")
	if(doUni)
		{
		unidirectionalFlag = sapply(1:nrow(data),FUN=function(x)
			{
			getMinStrand(data[x,,drop=FALSE],
				variantTypeCol=variantTypeCol,
				altCol=altCol,subName=subName)
			})
		} else {
		unidirectionalFlag = rep(NA,nrow(data))
		}
	data = cbind(data,unidirectionalFlag)
	uniCol = "unidirectionalFlag"
	# get germline filter
	print("set germline filter")
	if(doGermline)
		{
		if(germlineMeth=="sanger")
			{
			germlineFlag = sapply(1:nrow(data),
				FUN=function(x) getGermlineSanger(data[x,,drop=FALSE],
						variantTypeCol=variantTypeCol,
						altCol=altCol,
						refCol=refCol,
						subName=subName))
			} else {
			if(!is.na(germlineRefCountCol)&!is.na(germlineAltCountCol))
				{
				denom = (as.numeric(data[,germlineAltCountCol])+as.numeric(data[,germlineRefCountCol]))
				germlineFlag = as.numeric(data[,germlineAltCountCol])/denom
				} else {
				germlineFlag = rep(NA,times=nrow(data))
				}
			}
		} else {
		germlineFlag = rep(NA,times=nrow(data))
		}
	data = cbind(data,germlineFlag)
	germCol = "germlineFlag"
	# set patient
	print("set patient if missing")
	if(!patientCol%in%colnames(data))
		{
		patient = rep("patient",nrow(data))
		} else {
		patient = data[,patientCol]
		}
	# count number of patients with each variant
	print("count patients")
	variants = paste(data[,geneCol],data[,variantCol],sep=":")
	varTable = table(variants,patient)
	varTable[which(varTable>0)] = 1
	varTable = rowSums(varTable)
	patientCounts = varTable[variants]
	# return data
	return(list(data=data,civic=civic,sanger=sanger,
		mskcc=mskcc,blacklist=blacklist,
		genie=genie,uniCol=uniCol,germCol=germCol,
		patient=patient,
		civicMatch=civicMatch,sangerMatch=sangerMatch,
		mskccMatch=mskccMatch,blacklistMatch=blacklistMatch,
		genieMatch=genieMatch,patientCounts=patientCounts))
	}


# function to count the number of patients with the variant before prioritising
prioritiseVariant = function(
			tcgaCount=NULL,
			exacCount,
			impact,
			cadd,
			clinsig,
			variantClass,
			unidirectional,
			germline,
			civic,
			genie,
			sanger,
			mskcc,
			blacklist,
			doUni,doGermline,doCADD,
			patientCount
			)
	{
	INDEX <<- INDEX+1
	print(paste0(INDEX,". ",variantClass,":",patientCount))
	# classify variant priority
	priority = classifyVariant(
		  nPatients=patientCount,
                  tcgaCount=tcgaCount,
                  exacCount=exacCount,
		  impact=impact,
		  caddScore=cadd,
		  clinsig=clinsig,
		  variantClass=variantClass,
		  genie=genie,
		  sanger=sanger,
		  civic=civic,
		  mskcc=mskcc,
		  blacklist=blacklist,
		  unidirectionalFilter=unidirectional,
		  germlineFilter=germline,
		  doUni=doUni,doGermline=doGermline,doCADD=doCADD)
	return(c(priority,patientCount))
	}



# function to run priorities in full (pipeline)
runPriorities  = function(dataFile,genieFile,civicFile,sangerFile,mskccFile,blacklistFile,
		chromCol,posStartCol,posEndCol,tcgaCountCol,exacCountCol,
		geneCol,variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni=TRUE,doGermline=FALSE,doCADD=TRUE,
		germlineMeth="sanger",germlineRefCountCol=NA,germlineAltCountCol=NA,
		toRun="all",doWrite=TRUE,doParallel=TRUE)
	{
	# read and check data
	data = setupPrioritise(dataFile,genieFile,civicFile,sangerFile,mskccFile,blacklistFile,
		chromCol,posStartCol,posEndCol,tcgaCountCol,exacCountCol,
		geneCol,variantCol,impactCol,caddCol,clinsigCol,variantClassCol,
		patientCol,refCol,altCol,variantTypeCol,subName,
		doUni,doGermline,doCADD,germlineMeth,germlineRefCountCol,germlineAltCountCol)
	DATA <<- data
	# which variants to run
	if(toRun=="all") toRun = 1:nrow(data$data)
	# perform classification
	print("run prioritisation")
	INDEX<<-1
	if(!doParallel)
		{
		Priority = mapply(FUN=prioritiseVariant,
		    tcgaCount=data$data[toRun,tcgaCountCol],
		    exacCount=data$data[toRun,exacCountCol],
		    impact=data$data[toRun,impactCol],
		    cadd=data$data[toRun,caddCol],
		    clinsig=data$data[toRun,clinsigCol],
		    variantClass=data$data[toRun,variantClassCol],
		    unidirectional=data$data[toRun,data$uniCol],
		    germline=data$data[toRun,data$germCol],
		    civic=data$civicMatch,
		    sanger=data$sangerMatch,
		    genie=data$genieMatch,
		    mskcc=data$mskccMatch,
		    blacklist=data$blacklistMatch,
		    patientCount=data$patientCounts,
		    MoreArgs=list(doUni=doUni,doGermline=doGermline,doCADD=doCADD)
                    )
		} else {
		require(parallel)
		Priority = mcmapply(FUN=prioritiseVariant,
		    tcgaCount=data$data[toRun,tcgaCountCol],
		    exacCount=data$data[toRun,exacCountCol],
		    impact=data$data[toRun,impactCol],
		    cadd=data$data[toRun,caddCol],
		    clinsig=data$data[toRun,clinsigCol],
		    variantClass=data$data[toRun,variantClassCol],
		    unidirectional=data$data[toRun,data$uniCol],
		    germline=data$data[toRun,data$germCol],
		    civic=data$civicMatch,
		    sanger=data$sangerMatch,
		    genie=data$genieMatch,
		    mskcc=data$mskccMatch,
		    blacklist=data$blacklistMatch,
		    patientCount=data$patientCounts,
		    MoreArgs=list(doUni=doUni,doGermline=doGermline,doCADD=doCADD),
		    mc.cores=detectCores()
                    )

		}
	rownames(Priority) = c("priority","reason","patientCount")
	priorityTable = table(Priority["priority",])
	reasonTable = table(Priority["reason",])
	# combine data and results
	data = cbind(data$data,t(Priority))
	# write results
	if(doWrite)
		{
		outFile = getOutFile(dataFile)
		print("write results")
		write.table(data,
			file=outFile,
			quote=FALSE,row.names=FALSE,sep="\t")
		}
	return(data)
	}

# collate arguments into a list
collateArgs = function(args)
	{
	print("set arguments")
	args = list(
		# data files
		dataFile = args[1],
		genieFile = args[2],
		civicFile = args[3],
		sangerFile = args[4],
		mskccFile = args[5],
		blacklistFile = args[6],
		# column names
		chromCol = args[7],
		posStartCol = args[8],
		posEndCol = args[9],
		tcgaCountCol = args[10],
		exacCountCol = args[11],
		geneCol = args[12],
		variantCol = args[13], 
		impactCol = args[14],
		caddCol = args[15],
		clinsigCol = args[16],
		variantClassCol = args[17],
		patientCol = args[18],
		refCol = args[19],
		altCol = args[20],
		variantTypeCol = args[21],
		# variable names
		subName = args[22],
		# running options
		doUni = args[23],
		doGermline = args[24],
		doCADD = args[25],
		germlineMeth = args[26],
		germlineRefCountCol = args[27],
		germlineAltCountCol = args[28],
		doParallel = as.logical(args[29])
		)
	return(args)
	}



# prioritise
prioritiseVars = function(args)
	{
	argsList = collateArgs(args)
	do.call(runPriorities,argsList)
	}


# actually run prioritisation
prioritiseVars(args)

