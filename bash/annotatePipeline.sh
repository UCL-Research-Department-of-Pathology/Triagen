#!/bin/bash
# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
softDir="~/flanagan_chondrosarc"
vcfFile="vcf.vcf"
TCGAdir="/path/to/TCGA"
dataDir="/path/to/vcfDir"
outDir="/path/to/output"
genieFile=$softDir/data/genie/known_somatic_sites.bed

# get command line arguments
while getopts "h?svtdog:" opt; do
    case "$opt" in
    h|\?)
        echo "
	-s base directory of pipeline repository
	-v name of vcf file to annotate
	-t directory of TCGA and ICGC data
	-d directory that vcf file is in
	-o directory to write results to 
	-g bed file with known hotspot mutations from genie
	-h help
	"
        exit 0
        ;;
    s)  softDir=$OPTARG # base directory of pipeline repository
        ;;
    v)  vcfFile=$OPTARG # name of vcf file to annotate
        ;;
    t)  TCGAdir=$OPTARG # directory of TCGA and ICGC data
        ;;
    d)  dataDir=$OPTARG # directory that vcf file is in
        ;;
    o)  outDir=$OPTARG # directory to write results to 
        ;;
    g)  genieFile=$OPTARG # bed file with known hotspot mutations from genie
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

# run pipeline
baseFile=${vcfFile%.vcf}
# custom vcf2maf
perl $softDir/perl/vcf2maf.pl --input-vcf $dataDir/$vcfFile --output-maf $dataDir/$baseFile.maf
# add cadd score to maf
# v = vepFile
vepFile="$baseFile.vep.vcf" 
python $softDir/py/add_cadd.py -i $dataDir/$baseFile.maf -v $dataDir/$vepFile -o $outDir/$baseFile-withCadd.maf
# get TCGA info
python $softDir/py/comparing_maf_files.py -i $outDir/$baseFile-withCadd.maf -a $TCGAdir -o $outDir/$baseFile-withCadd-withTCGA.maf
# prioritise variants based on filters
Rscript --vanilla $softDir/R/prioritiseVariants.R $outDir/$baseFile-withCadd-withTCGA.maf $genieFile
# remove intermediates
rm $dataDir/$baseFile.maf
rm $outDir/$baseFile-withCadd.maf
rm $outDir/$baseFile-withCadd-withTCGA.maf

