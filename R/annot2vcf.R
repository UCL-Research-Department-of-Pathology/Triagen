annot2vcf = function(data,outFile,chromCol=5,posCol=6,idCol=NULL,refCol=7,altCol=8,qualCol=9,filterCol=10,flankingCol="flanking_bps",split=1,idAddition="id")
	{
	if(is.null(idCol))
		{
		id = paste0(idAddition,1:nrow(data))
		data = cbind(data,id)
		idCol = ncol(data)
		}
	# create vcf
	vcf = data[,c(chromCol,posCol,idCol,refCol,altCol,qualCol,filterCol)]
	# deal with indels
	indexDel = grep("-",vcf[,5])
	if(length(indexDel)>0)
		{
		preceding = sapply(data[indexDel,flankingCol],FUN=function(x) strsplit(x,split="")[[1]][2])
		vcf[indexDel,4] = paste0(preceding,vcf[indexDel,4])
		vcf[indexDel,5] = preceding
		vcf[indexDel,2] = vcf[indexDel,2]-1
		}
	indexIns = grep("-",vcf[,4])
	if(length(indexIns)>0)
		{
		preceding = sapply(data[indexIns,flankingCol],FUN=function(x) strsplit(x,split="")[[1]][1])
		vcf[indexIns,5] = paste0(preceding,vcf[indexIns,5])
		vcf[indexIns,4] = preceding
		vcf[indexIns,2] = vcf[indexIns,2]-1
		}
	# column headings
	newHeads = c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER")
	vcf = rbind(newHeads,vcf)
	# create vcf headers
	filters = unique(unlist(sapply(data[,filterCol],FUN=function(x) strsplit(x,split=";")[[1]])))
	headers = c('##fileformat=VCFv4.0',
		paste0('##fileDate=',gsub("-","",Sys.Date())),
		paste0('##FILTER=<ID=',filters,',Description="',filters,'">'))
	headFillers = matrix("",ncol=ncol(vcf)-1,nrow=length(headers))
	headers = cbind(headers,headFillers)
	colnames(headers) = colnames(vcf)
	# combined headers and data
	if(split>1)
		{
		# get file names
		ending = rev(strsplit(outFile,split="[.]")[[1]])[1]
		fileBase = gsub(paste0("[.]",ending),"",outFile)
		outFiles = paste0(fileBase,"-split-",1:split,".",ending)
		# split data
		binSize = ceiling(nrow(vcf)/split)
		cohorts = 1:nrow(vcf)
		cohorts = split(cohorts, ceiling(seq_along(cohorts)/binSize))
		vcfs = lapply(cohorts,FUN=function(x) rbind(headers,vcf[x,,drop=FALSE]))
		# write vcfs
		sapply(1:length(vcfs),FUN=function(x) write.table(vcfs[[x]],sep="\t",file=outFiles[x],row.names=FALSE,col.names=FALSE,quote=FALSE))
		} else {
		vcf = rbind(headers,vcf)
		# write vcf
		write.table(vcf,sep="\t",file=outFile,row.names=FALSE,col.names=FALSE,quote=FALSE)
		}
	}


#data = read.table("MPNST_good_Pind_Cave.txt",as.is=TRUE)
#annot2vcf(data,"MPNST_good_Pind_Cave.vcf")
