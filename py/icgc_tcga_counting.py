import os, re, sys
import subprocess
import argparse

def read_anno(anno_file):
	anno = {}
	patients_tmp = {}
	patient_count = {}
	with open(anno_file) as f:
		for line in f:
			line= line.rstrip()
			word = line.split('\t')
			name = word[0]
			if name not in anno:
				anno[name] = {}
			if name not in patient_count:
				patient_count[name] = int(word[5])
			if (word[1], word[2], word[3], word[4]) not in anno[name]:
				anno[name][(word[1], word[2], word[3], word[4])] = int(word[5])
	return anno, patient_count


def read_input_file(maf, anno, patient_dict, outfile):
	output = open(outfile, 'w')
	with open(maf) as f:
		next(f)
		header = next(f)
		header = header.rstrip()
		output.write(header),
		head = header.split('\t')
		alt_col = head.index('Tumor_Seq_Allele2')
		chr_col = head.index('Chromosome')
		try:
			start_col = head.index('Start_position')
		except ValueError as e:
			start_col = head.index('Start_Position')
		try:
			end_col = head.index('End_position')
		except ValueError as j:
			end_col = head.index('End_Position')
		output.write('\tNo of Occurances in TCGA/ICGC'),
		for s in sorted(anno.keys()):
			output.write('\t{}'.format(s)),
		output.write('\n'),
		for line in f:
			count = {}
			line = line.rstrip(' \n')
			word = line.split('\t')
			alt_col = head.index('Tumor_Seq_Allele2')
			chr_col = head.index('Chromosome')
			try:
				start_col = head.index('Start_position')
			except ValueError as e:
				start_col = head.index('Start_Position')
			try:
				end_col = head.index('End_position')
			except ValueError as j:
				end_col = head.index('End_Position')
			total_count = 0
			for study in anno:
				if (word[chr_col], word[start_col], word[end_col], word[alt_col]) in anno[study]:
					count[study] = anno[study][(word[chr_col], word[start_col], word[end_col], word[alt_col])]
					total_count += anno[study][(word[chr_col], word[start_col], word[end_col], word[alt_col])]
					count[study] = float(count[study])/patient_dict[study]
			output.write(line),
			output.write('\t{}'.format(total_count)),
			for c in sorted(anno.keys()):
				if c in count:
					output.write('\t{}'.format(count[c])),
				else:
					output.write('\t'),
			output.write('\n'),
	output.close()

def main():
	parser = argparse.ArgumentParser(description='Annotates a maf file to another\n')
	parser.add_argument('-i','--input', help='Input maf file', required=False)
	parser.add_argument('-a','--anno', help='TCGA/ICGC summary file', required=False)
	parser.add_argument('-o','--output', help='Output file', required=False)
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = vars(parser.parse_args())

	anno, patient_count = read_anno(args['anno'])
	read_input_file(args['input'], anno, patient_count, args['output'])

main()
