import os
import vcf
import shutil
import sys
import re

if os.name == 'nt':
	GITHUB_FOLDER = os.path.join(os.getenv('USERPROFILE'), 'Documents', 'Github')
else:
	GITHUB_FOLDER = os.path.join(os.getenv('HOME'), 'Documents', 'Github')
sys.path.append(GITHUB_FOLDER)

import pytools.systemtools as systemtools
import pytools.filetools as filetools
import varianttools.callertools as callertools
from pprint import pprint


def copyVcf(source, destination):
	with open(source, 'r') as input_file:
		reader = vcf.Reader(input_file)
		if 'Varscan' in source:
			reader.formats['DP4'] = reader.formats['DP4']._replace(num=4)
			reader.formats['DP4'] = reader.formats['DP4']._replace(type='Integer')
		with open(destination, 'w') as output_file:
			writer = vcf.Writer(output_file, reader)
			for record in reader:
				filterOut = '/' in str(record.ALT[0]) or '/' in record.REF
				if not filterOut:
					try:
						writer.write_record(record)
					except ValueError:
						print(record)
	return destination


def splitVcfByChromosome(source, output_folder, create_subfolders = False):
	""" Separates a vcf file into separate files for each chromosome.
		Assumes the file is sorted.
		
		Parameters
		----------
			source: string [PATH]
			output_folder: string [PATH]
			create_subfolders: bool; default False
				If 'True', each chromosome will be saved to a separate folder.
	"""
	basename = os.path.basename(source)
	basename, ext = os.path.splitext(basename)
	_match_chroms = "chr[0-9MT]{1,3}$"
	_match_chroms = re.compile(_match_chroms)

	with open(source, 'r') as input_vcf_file:
		reader = vcf.Reader(input_vcf_file)
		#pprint(reader.contigs)
		chromosomes = {i:list() for i in reader.contigs}
		# Sort the records by chromosome
		for record in reader:
			chrom = record.CHROM
			if chrom not in chromosomes: 
				chromosomes[chrom] = list() 
			chromosomes[chrom].append(record)

		for chromosome, record_list in chromosomes.items():
			match = _match_chroms.search(chromosome)
			if not match: continue
			print(chromosome, match)
			output_basename = "{}.{}.vcf".format(basename, chromosome)
			print(output_basename)
			if create_subfolders:
				chromosome_folder = os.path.join(output_folder, chromosome)
			else: 
				chromosome_folder = output_folder

			output_filename       = os.path.join(chromosome_folder, output_basename)

			filetools.checkDir(chromosome_folder, True)
			with open(output_filename, 'w') as output_vcf:
				writer = vcf.Writer(output_vcf, reader)
				if len(record_list) > 0:
					for record in record_list:
						writer.write_record(record)


def splitCallsetByChromosome(callset, output_folder):
	""" Splits the files in a patients callset.

		Parameters
		----------
			callset: string, dict<>
				Eith a folder containing a callet, or
				a pre-processed callset.
			output_folder: string
				The folder to save the split files to.
				Each chromosome will be saved to a separate sub-folder.
	"""
	filetools.checkDir(output_folder)

	if isinstance(callset, str):
		classifier = callertools.CallerClassifier()
		callset = classifier(callset, key = 'original')

	for caller, filename in callset.items():
		splitVcfByChromosome(filename, output_folder, create_subfolders = True)


def splitVcf(filename, output_folder = None):
	""" Separates a vcf file into indel and snp sections.
		Parameters
		----------
			filename: string [PATH]
				THe vcf file to split.
			output_folder: path
				The folder to save the output files to. If None,
				the files will be saved to the folder the original file
				resides.

		Returns
		-------
			output: dict<>
				* 'indel': path
					The indel file
				* 'snp': path
					The snp file
	"""

	path, basename = os.path.split(filename)
	if output_folder is None: output_folder = path
	basename, ext = os.path.splitext(basename)
	snp_filename = os.path.join(output_folder, basename + ".snp.vcf")
	indel_filename=os.path.join(output_folder, basename + ".indel.vcf")
	
	with open(filename, 'r') as vcf_file:
		reader          = vcf.Reader(vcf_file)

		snp_writer      = vcf.Writer(open(snp_filename,   'w'), reader)
		indel_writer    = vcf.Writer(open(indel_filename, 'w'), reader)
		for record in reader:
			if record.is_snp:
				snp_writer.write_record(record)
			elif record.is_indel:
				indel_writer.write_record(record)
			else:
				print("vcftools.splitVcf: ", record)

	result = {
		'snp': snp_filename,
		'indel': indel_filename

	}
	return result


def splitCallset(callset, output_folder, **kwargs):
	split_callset = dict()
	#pprint(callset)
	for caller_name, source in callset.items():
		#print("\tsplitting ", caller_name)
		#print("\t\tSource: ", source)
		
		if 'indel' in caller_name or 'snp' in caller_name:
			destination = os.path.join(output_folder, os.path.basename(source))
			shutil.copy2(source, destination)
			split_callset[caller_name] = destination
		else:
			result = splitVcf(source, output_folder)
			split_callset[caller_name + '-indel'] = result['indel']
			split_callset[caller_name + '-snp']   = result['snp']
			#print("\t\tDestination: ", result['indel'])
			#print("\t\tDestination: ", result['snp'])

	return split_callset


def fixCallerOutputs(callset, somaticseq_folder, **kwargs):
	"""
		Required Parameters
		-------------------
			variants: dict<str:str> [dict<caller_name: caller_output>]
				The callset to fix
			somaticseq_folder: str [path]
				Th folder containing the somaticseq program.
		Optional Parameters
		-------------------
			sample: dict<>
			patientId: str

	"""

	modify_vjsd_script   = os.path.join(somaticseq_folder, "modify_VJSD.py")
	
	fixed_callset = dict()
	for caller, source in callset.items():
		if 'output_folder' in kwargs:
			output_folder = kwargs['output_folder']
		else:
			output_folder = os.path.dirname(source)

		if 'patientId' in kwargs:
			basename = "{}.{}.corrected.vcf".format(kwargs['patientId'], caller)
		else:
			basename = os.path.basename(source)
			basename, ext = os.path.splitext(basename)
			basename = "{}.corrected.vcf".format(basename)

		destination = os.path.join(output_folder, basename)
		if 'varscan' in caller:
			command = """python3 {program} -method VarScan2 -infile {infile} -outfile {outfile}"""
		elif 'somaticsniper' in caller:
			command = """python3 {program} -method SomaticSniper -infile {infile} -outfile {outfile}"""
		elif 'muse' in caller:
			command = """python3 {program} -method MuSE -infile {infile} -outfile {outfile}"""
		else:
			command = None

		if command:
			command = command.format(
				program = modify_vjsd_script,
				infile  = source,
				outfile = destination)
			systemtools.Terminal(command, use_system = True)
		else:
			shutil.copy2(source, destination)
		fixed_callset[caller] = destination

	return fixed_callset

if __name__ == "__main__":
	source_folder = "/home/upmc/Documents/Data/raw_snp_output/TCGA-2H-A9GR/"
	output_folder = "/home/upmc/Documents/Data/raw_snp_output/TCGA-2H-A9GR-CHROMS/"
	splitCallsetByChromosome(source_folder, output_folder)


