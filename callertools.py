""" Tools for organizing the output from each caller """
import os
import re
import sys
from pprint import pprint
if os.name == 'nt':
	GITHUB_FOLDER = os.path.join(os.getenv('USERPROFILE'), 'Documents', 'Github')
else:
	GITHUB_FOLDER = os.path.join(os.getenv('HOME'), 'Documents', 'Github')
sys.path.append(GITHUB_FOLDER)

import pytools.filetools as filetools
import pytools.tabletools as tabletools
import progressbar


class CallerClassifier:
	""" Classifies a callset based on caller via regexes.
	"""
	def __init__(self):
		self.regexes = self._defineRegexes()

	def __call__(self, folder, **kwargs):
		"""
			Keyword Arguments
			-----------------
				'type': {'indel', 'snp'}
		"""
		all_files = filetools.listAllFiles(folder, **kwargs)
		#from pprint import pprint
		#pprint(all_files)
		results = dict()
		for filename in all_files:
			caller = self._classifyFilename(filename.lower())
			if caller:
				if 'indel' in caller or 'snp' in caller:
					_addition = ''
				elif 'indel' in filename:
					_addition = '-indel'
				elif 'snp' in filename or 'snv' in filename:
					_addition = '-snp'
				else:
					_addition = ''
				caller_key = caller[:-1] + _addition
				results[caller_key] = filename

		if 'type' in kwargs:
			results = {k: v for k, v in results.items() if kwargs['type'] in k}

		return results

	def _classifyFilename(self, filename):
		filename = filename.lower()
		for caller_base_name, templates in self.regexes.items():
			for caller, regex in templates.items():
				match = regex.search(filename)
				if match:
					return caller     
		

	@staticmethod
	def _defineRegexes():
		muse_regex = {
			'muse0': ".muse\..*\.?vcf$",
		}
		mutect_regex = {
			'mutect20': "\.mutect2\..*\.?vcf$"
		}
		somaticsniper_regex = {
			'somaticsniper0': "\.somaticsniper.*\.vcf$"
		}
		strelka_regex = {
			'strelka-indel0': "\.passed\.somatic\.indels\..*\.?vcf$",
			'strelka-snp0': "\.passed\.somatic\.sn[vp]s\..*\.?vcf$"
		}
		varscan_regex = {
			'varscan-indel0': "\.raw\.indel\..*\.?vcf$",
			'varscan-snp0':   "\.raw\.snp\.somatic\.hc\..*\.?vcf$",
			'varscan-snp1': ".varscan.snp.somatic.hc.vcf$"
		}
		haplotypecaller_regex = {
			'haplotypecaller-rna0': "haplotypecaller.rna.filtered_variants.vcf$",
			'haplotypecaller-dna0': "haplotypecaller.rna.raw_variants.vcf$"
		}
		# TCGA-2H-A9GF-CHR1-11A_vs_TCGA-2H-A9GF-CHR1-01A.varscan.snp.Somatic.hc
		regexes = {
			'haplotypecaller': 	{k:re.compile(v) for k, v in haplotypecaller_regex.items()},
			'muse': 			{k:re.compile(v) for k, v in muse_regex.items()},
			'mutect2': 			{k:re.compile(v) for k, v in mutect_regex.items()},
			'somaticsniper': 	{k:re.compile(v) for k, v in somaticsniper_regex.items()},
			'strelka': 			{k:re.compile(v) for k, v in strelka_regex.items()},
			'varscan': 			{k:re.compile(v) for k, v in varscan_regex.items()}
		}
		return regexes


def getCallerStatus(folder):
		"""
			Parameters
			----------
				folder: Folder containing all patient_folders.
		"""
		classifier = CallerClassifier()
		expected_callers = ['muse', 'mutect2', 'somaticsniper', 
			'strelka-indel', 'strelka-snp', 'varscan-indel', 'varscan-snp']
		patient_barcodes = sorted(os.listdir(folder))
		results = list()
		print("Parsing folder...")
		bar = progressbar.ProgressBar(max_value = len(patient_barcodes))

		for index, patient_barcode in enumerate(patient_barcodes):
			bar.update(index)
			abs_path = os.path.join(folder, patient_barcode)
			callset = classifier(abs_path)
			result = {
				'patientId': patient_barcode,
			}
			for caller_name in expected_callers:
				result[caller_name] = caller_name in callset

			all_patient_files = filetools.listAllFiles(abs_path)
			total_size = sum([os.path.getsize(fn) for fn in all_patient_files])
			result['totalSize'] = total_size
			results.append(result)

		results = tabletools.Table.fromList(results)

		return results


class GATKMergeSampleCallsets:
	"""Merges callsets via GATK CombineVariants.
		WARNING: CombineVariants._copy_vcf() uses hard filter to remove variants with '/' formatting
		Usage
		-----
			merger = MergeSampleCallsets()
			output = merger(sample, output_folder, callset)
		Requirements
		------------
			GATK
			Reference Genome
	"""
	def __init__(self, merge_options = None, **kwargs):
		""" Parameters
			----------
				sample:
				merge_options: the general options being used for the genomics pipeline.
				output_folder:
				variants:
		"""
		#log_message = "GATKMergeSampleCallsets(merge_options = {})".format(merge_options)
		#print(log_message)
		#print('\tKeyword Arguements:')
		#for k, v in kwargs.items():
		#	print("\t{}\t{}".format(k, v))

		if merge_options is not None:
			self.gatk_program = merge_options['Programs']['GATK']
			self.reference = merge_options['Reference Files']['reference genome']
		else:
			self.gatk_program = kwargs.get('GATK', kwargs.get('program'))
			self.reference = kwargs.get('reference')

		if self.gatk_program is None or self.reference is None:
			message = "Missing a dependancy: "
			if self.gatk_program is None: message += ', GATK'
			if self.reference is None: message += ', reference'
			raise ValueError(message)

	def __call__(self, callset, **kwargs):
		""" Merges the provided callset
			Keyword Arguments
			-----------------
				* 'filename': The filename to save the merged callset as.
				* 'kind': If filename is not provided, 'kind' indicates the patient
					folder to use. Both 'patientId' and 'kind' must be provided.

		"""
		output_filename = kwargs.get('filename')
		self._checkCallsetFormat(callset)
		output_filename = self.gatkCombineVariants(callset, output_filename)
		return output_filename
			
	def _modify_merged_vcf(self, filename):
		""" Adds the VAF to the 'Info' field of the output file.
			The new file will be saved to the same folder as the original.
		:param filename: string
			Path to the merged file.
		:return: string
			path to the output file.
		"""
		output_folder = os.path.dirname(filename)
		basename = os.path.splitext(os.path.basename(filename))[0]
		basename = basename + ".modified.vcf"
		output_file = os.path.join(output_folder, basename)

		with open(filename, 'r') as vcf_file:
			reader = vcf.Reader(vcf_file)
			reader.infos['VAF'] = reader.formats['FREQ']._replace(type='Float')
			
			with open(output_file, 'w') as file2:
				writer = vcf.Writer(file2, reader)
				for record in reader:
					VAF = self._getVAF(record)
					record.INFO['VAF'] = VAF
					writer.write_record(record)

		return output_file

	@staticmethod
	def _checkCallsetFormat(callset):
		""" Ensures a callset is formatted as a dictionary with a single file
			per caller. 
		"""
		caller_name_format_status = any('-' in c for c in callset.keys())

		callset_valid = not caller_name_format_status
		if not callset_valid:
			message = "The callset names include '-'"
			raise ValueError(message)
		return caller_name_format_status

	@staticmethod
	def _getVAF(record):
		""" Use DP4 instead of DP
			Parameters
			----------
				record: 

			Returns
			-------
				result: dict<>
					* 'alleles': The number of alternate alleles.
					* 'reads': The number of reads used to calculate the vaf.
					* 'vaf': The variant allele frequency, in the range 0 - 100.
		"""
		info_keys = record.INFO.keys()
		record_sample = [i for i in record.samples if i.sample == 'TUMOR'].pop()
		if 'FREQ' in info_keys:
			VAF = float(record_sample['FREQ'][:-1])
		elif 'DP' in info_keys and 'AD' in info_keys:
			reads = record_sample['DP']
			alleles = record_sample['AD'][1]
			VAF = (alleles / reads) * 100
		elif 'AF' in info_keys:
			VAF = record_sample['AF']
		elif 'DP4' in info_keys:
			reads = sum(record_sample['DP4'])
			alleles = sum(record_sample['DP4'][2:])
			VAF = alleles / reads
		else:
			alleles = [i for i in ['A', 'C', 'G', 'T'] if i != record.REF]
			sample_reads = sum([record_sample[i+'U'][1] for i in (alleles + [record.REF])])
			sample_alleles = sum([record_sample[i+'U'][1] for i in alleles])
			if sample_reads == 0: sample_vaf = 0
			else:
				sample_vaf = sample_alleles / sample_reads
			VAF = sample_vaf
		return VAF

	def _modify_variants(self, vcfs):
		""" Some caller outputs are inconsistent and need to be modified.
		"""
		output_variants = dict()
		for caller, vcf_filename in vcfs.items():
			current_output_folder = os.path.dirname(vcf_filename)
			basename = os.path.splitext(os.path.basename(vcf_filename))[0] + ".modified.vcf"
			output_file = os.path.join(current_output_folder, basename)

			with open(vcf_filename, 'r') as vcf_file:
				reader = vcf.Reader(vcf_file)
				if 'varscan' in caller:
					reader = self._modify_varscan_output(reader)
				output_variants[caller] = self._copy_vcf(reader, output_file)

		return output_variants
	@staticmethod
	def _copy_vcf(reader, output_file):
		with open(output_file, 'w') as file1:
			writer = vcf.Writer(file1, reader)
			for record in reader:
				filterOut = '/' in str(record.ALT[0]) or '/' in record.REF
				if filterOut:
					pass
				else:
					writer.write_record(record)
		return output_file

	@staticmethod
	def _modify_varscan_output(reader):
		reader.formats['DP4'] = reader.formats['DP4']._replace(num=4)
		reader.formats['DP4'] = reader.formats['DP4']._replace(type='Integer')
		return reader

	def _combineSplitVariants(self, patientId, output_folder, callset):
		""" Merges the indel/snp variants separately.
			Parameters
			----------
				sample: dict<>, string
					sample information as a dict or the patient id.
				output_folder: path
					Folder to save the merged variant files in.
				callset: Callset; default None
					If not provided, a generic callset will be created.
		"""
		if isinstance(patientId, dict):
			patientId = patientId['PatientID']

		snp_variants   = callset('snp')
		indel_variants = callset('indel')

		basename = "{}.merged_variants".format(patientId)
		output_snp_filename   = os.path.join(output_folder, basename + '.snp.vcf')
		output_indel_filename = os.path.join(output_folder, basename + '.indel.vcf')

		merged_snp_filname    = self.gatkCombineVariants(snp_variants,   output_snp_filename)
		merged_indel_filename = self.gatkCombineVariants(indel_variants, output_indel_filename)

		result = {
			'snp': merged_snp_filname,
			'indel': merged_indel_filename
		}
		return result

	def gatkCombineVariants(self, variants, output_file):
		""" Uses GATK CombineVariants to merge the calls from each caller into a single file.
			Parameters
			----------
				variants: dict<caller, path>
					A dictionary linkng each caller to its harmonized output.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.harmonized.vcf
			Returns
			-------
				Output_file: string
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.merged.vcf
		"""

		order = "mutect2,varscan,strelka,muse,somaticsniper"  # ordered by VAF confidence
		variant_command = ['--variant:{} "{}"'.format(k, v) for k, v in variants.items()]
		variant_command = ' \\\n'.join(variant_command)
		command = """java -jar "{gatk}" \
			-T CombineVariants \
			-R "{reference}" \
			{variants} \
			-o "{output}" \
			-genotypeMergeOptions PRIORITIZE \
			-priority {rod}"""
		command = command.format(
				gatk =      self.gatk_program,
				reference = self.reference,
				variants =  variant_command,
				rod =       order,
				output =    output_file)

		systemtools.Terminal(command)
		return output_file

	def catVariants(self, left, right):
		""" Combines the SNV and Indel files. Assumes both are saved in the same folder. """
		l = os.path.splitext(os.path.splitext(left)[0])[0]
		output_file = l + '.cat.vcf'

		command = """java -cp {GATK} org.broadinstitute.gatk.tools.CatVariants \
			-R {reference}\
			-V {left} \
			-V {right} \
			-out {output}
		""".format(
			GATK = self.gatk_program,
			reference = self.reference,
			left = left,
			right = right,
			output = output_file)
		systemtools.Terminal(command)

		return output_file

def classify(folder):
	""" Retrieves the callset for each caller """
	return _CLASSIFIER(folder)

def merge(callset, filename, program, reference):
	"""
		Merges a callset.
		Parameters
		----------
			callset
			filename
			program: Path to GATK
			reference: path to the reference genome.
	"""
	gatk = GATKMergeSampleCallsets(GATK = program, referece = reference)
	result = gatk(
		callset = callset,
		filename = filename,
	)
	return result


_CLASSIFIER = CallerClassifier()

if __name__ == "__main__":
	#test_folder = "G:\\Data\\TCGA-ESCA\\raw_snp_output\\"
	#test_folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/3_called_variants/"
	test_folder = "/media/upmc/LMD_boot/home/upmc/Documents/Variant_Discovery_Pipeline/3_called_variants"
	table = getCallerStatus(test_folder)
	print(table.df)
	table.df.to_csv("LMD_caller_status.tsv", sep = "\t")
	print("Saved!")
