# Triagen
Strategy for triaging genomic variants

For full details of the method and citation, please see [Steele et al. (2019)](https://www.sciencedirect.com/science/article/pii/S1535610819300972).

Depends on: VEP (https://github.com/Ensembl/ensembl-vep)

Attempt to automatically triage variants into those that we think are importnat, and those that we think are artefact/unimportant.

General workflow:
	1 - Annotate variants with VEP

	2 - Count occurences of variants in TCGA/ICGC

	3 - subset to regions of interest

	4 - remove variants that fail quality filters

	5 - filter variants believed to be artefact (germline/unidirectional)

	6 - filter variants that are common in population (ExAC)

	7 - cross reference with previous datasets (Genie,CIVIC,Sanger)

	8 - cross reference with TCGA

	9 - refine importance with predicted pathogenicity (CADD)

./py/run_new_anno performs steps 1 and 2.

./R/filterNonCoding.R performs steps 3 and 4.

./R/prioritiseVariants performs steps 5-9.
