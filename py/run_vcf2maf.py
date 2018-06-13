import os, subprocess, re, sys
import argparse


def anno(ifile, vep_path, outfile):
	commands =  '/home/regmplo/Scratch/18_chondrosarcoma/flanagan_chondrosarc/perl/vcf2maf.pl  --vep-path {0} --input-vcf {1}  --output-maf {2}.maf '\
		' --vep-data {0}/data --ref-fasta {0}/data/homo_sapiens/87_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa  '\
		' \n'.format(vep_path, ifile, outfile )
	qsub(commands, 'tmp', 'anno', 1, None, '2G', '/home/regmplo/Scratch/18_chondrosarcoma/undiff_annotation/logs/')
	#subprocess.call(commands.split())

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
#$ -l h_rt=48:0:0
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
source ~/env3/bin/activate
module unload perl/5.22.0
module load perl/5.16.0
{0}'''.format(command)
	return wrap


def main():
	parser = argparse.ArgumentParser(description='Annotates a maf file to another\n')
	parser.add_argument('-i','--input', help='Input vcf', required=False)
	parser.add_argument('-v','--vep_path', help='Vep path to installation', required=False)
	parser.add_argument('-o','--output', help='Output file', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	anno(args['input'], args['vep_path'],  args['output'])

main()