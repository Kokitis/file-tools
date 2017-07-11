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
		# TCGA-2H-A9GF-CHR1-11A_vs_TCGA-2H-A9GF-CHR1-01A.varscan.snp.Somatic.hc
		regexes = {
			'muse': {k:re.compile(v) for k, v in muse_regex.items()},
			'mutect2': {k:re.compile(v) for k, v in mutect_regex.items()},
			'somaticsniper': {k:re.compile(v) for k,v in somaticsniper_regex.items()},
			'strelka': {k:re.compile(v) for k, v in strelka_regex.items()},
			'varscan': {k:re.compile(v) for k, v in varscan_regex.items()}
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


if __name__ == "__main__":
	#test_folder = "G:\\Data\\TCGA-ESCA\\raw_snp_output\\"
	#test_folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/3_called_variants/"
	test_folder = "/media/upmc/LMD_boot/home/upmc/Documents/Variant_Discovery_Pipeline/3_called_variants"
	table = getCallerStatus(test_folder)
	print(table.df)
	table.df.to_csv("LMD_caller_status.tsv", sep = "\t")
	print("Saved!")
