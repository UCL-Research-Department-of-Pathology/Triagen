#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test arguments: if not, return an error
if (!length(args)%in%c(1,4) 
	{
	stop("Incorrect number of arguments", call.=FALSE)
	} else if (length(args)!=1) {
	# default column headings
	args[2:4] = c(
		"Variant_Classification",	# Effect of variant
		"CLPM",				# median clipped reads
		"ASMD"				# median alignment score
		)
	}


# column names for columns of interest
effectCol = args[2]
CLPMcol = args[3]
ASMDcol = args[4]
# load data
dataFile = args[1]
data = read.table(dataFile,head=TRUE,sep="\t",comment.char="@",quote="",as.is=TRUE)

# remove intronic etc
#removeIndex = which(data[,effectCol] %in% c(
#	"3'Flank","3'UTR","5'Flank","5'UTR",
#	"IGR","Intron","RNA"))
print(dim(data))
removeIndex = which(data[,effectCol] %in% c(
	"IGR","Intron","RNA"))
data = data[-removeIndex,]
print(dim(data))

# remove low quality calls
ASMD = as.numeric(data[,ASMDcol])
keepIndex = which(is.na(ASMD)|ASMD>=140)
data = data[keepIndex,]
print(dim(data))
CLPM = as.numeric(data[,CLPMcol])
keepIndex = which(is.na(CLPM)|CLPM==0)
data = data[keepIndex,]
print(dim(data))

# get out file name
splitName = strsplit(dataFile,split="[.]")[[1]]
fileEnding = rev(splitName)[1]
fileName = splitName[1]
outFile = paste0(fileName,"-filteredNonCoding.",fileEnding)
# write results
if(fileEnding=="csv")
	{
	write.csv(data,
		file=outFile,
		quote=FALSE,row.names=FALSE,na="")
	} else {
	write.table(data,
		file=outFile,
		quote=FALSE,row.names=FALSE,sep="\t",na="")
	}









