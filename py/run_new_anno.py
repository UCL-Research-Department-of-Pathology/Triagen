
import os, subprocess, re, sys

pipelineDir = '/home/ucbtcds/pipelines/Triagen'
vep_path = '/home/ucbtcds/software/ensembl-vep'
vepcache_path = '/home/ucbtcds/Scratch/vepCache'
exacFile = '/home/ucbtcds/Scratch/data/ExAC/ExAC.r0.3.1.sites.vep.vcf.gz'
logDir = '/home/ucbtcds/Scratch/logs/11_MPNST/'
tcgaIcgcFile = '/home/ucbtcds/Scratch/data/TCGA/exomes/icgc_tcga_data_summary.txt'

def anno(idir):
	ifiles = [f for f in os.listdir(idir) if f.endswith('.vcf')]
	for i in ifiles:
		commands =   '{4}/perl/custom_vcf2maf.pl --vep-path {0} --input-vcf {1}/{2}  --output-maf {1}/{2}.maf '\
			' --vep-data {3} --ref-fasta {3}/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz '\
			' --tumor-id TUMOR --normal-id NORMAL --cache-version 90 '\
			' --filter-vcf {5}'.format(vep_path, idir, i, vepcache_path, pipelineDir, exacFile )
		#print commands
		qsub(commands, 'anno', i, 1, None, '4G', logDir)
		#subprocess.call(commands, shell=True)

def anno2(idir):
	ic_sc = pipelineDir + '/py/icgc_tcga_counting.py'
	a_file = tcgaIcgcFile
	indir = os.getcwd()
	ifiles = [f for f in os.listdir(idir) if f.endswith('.vcf.maf')]
	for i in ifiles:
		command =  'python {0} -i {1}/{3} -a {2} -o {1}/{3}.icgc.maf\n'.format(ic_sc, idir, a_file, i)
		qsub(command, 'anno2', i, 1, None, '1G', logDir)

def qsub(command, sample, name, threads, space, mem, outdir):
	wrapper = submit(command, sample, name, threads, space, mem, outdir)
	output = open('{}/{}_{}.sh'.format(outdir, sample, name), 'w')
	output.write(wrapper)
	output.close()
	stdout = subprocess.Popen('qsub {}/{}_{}.sh'.format(outdir, sample, name), shell=True, stdout=subprocess.PIPE).communicate()[0]
	stdout = re.sub('Your job', '', stdout)
	stdout2 = stdout.split()
	print stdout2[0]

def submit(command, sample, name, threads, space, mem, outdir):
	wrap = '''#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=3:50:0
#$ -l mem={}
'''.format(mem)
	if space:
		wrap = wrap + '#$ -l tmpfs={}'.format(space)
	wrap = wrap + '''
#$ -N {0}_{1}
#$ -wd {2}
'''.format(sample, name, outdir, mem)
	if int(threads) > 1:
		wrap = wrap + '#$ -pe mpi {}'.format(threads)
	#wrap = wrap + '#$ -pe smp {}'.format(threads)
	wrap = wrap + '''
module load samtools
module load perlbrew
source ~/.bashrc
{0}'''.format(command)
	return wrap


def main():
	anno2('/home/ucbtcds/Scratch/data/11_MPNST') #Annotate test calls


main()
