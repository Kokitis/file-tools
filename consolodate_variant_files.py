import os
import re
import shutil
""" Transers the output from each caller to a new folder. Gathers statistics on which callers have successfully been run for each case """

base_folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/" #Base folder of the variant discovery pipeline

def search_folder(folder, string):
	matches = [fn for fn in os.listdir(folder) if string in fn]
	if len(matches) == 0:
		result = None
	elif len(matches) > 1:
		#print("Found mutliple matching values: ", matches)
		result = matches
	else:
		result = matches
	if result is not None:
		result = [os.path.join(folder, fn) for fn in result]
	return result

def get_folders(source_folder, patient_barcode):
	muse_folder = os.path.join(source_folder, patient_barcode, 'MuSE')
	mutect_folder = os.path.join(source_folder, patient_barcode, 'MuTect2')
	somaticsniper_folder = os.path.join(source_folder, patient_barcode, 'SomaticSniper')
	strelka_folder = os.path.join(source_folder, patient_barcode, 'Strelka', 'results')
	varscan_folder = os.path.join(source_folder, patient_barcode, 'Varscan')
	return muse_folder, mutect_folder, somaticsniper_folder, strelka_folder, varscan_folder

def transferVariants(source_folder, destination_folder, DEBUG = False):
	""" Transfers usable somatic vcf files from one folder to another.
		Assumes the directory is structured as:
			<patient barcode>/<caller name>/<variant files>
		Parameters
		----------
			source_folder: string [Folder Path]
				A folder containing variants for patients.
			destination_folder: string [Folder Path]
				The folder to move the files to.
	"""

	patient_barcodes = [fn for fn in os.listdir(source_folder) if '-' in fn]

	for index, patient_barcode in enumerate(patient_barcodes):
		print("{0} of {1}".format(index+1, len(patient_barcodes)))
		muse_folder, mutect_folder, somaticsniper_folder, strelka_folder, varscan_folder = get_folders(source_folder, patient_barcode)
		#muse_folder = os.path.join(source_folder, patient_barcode, 'MuSE')
		#mutect_folder = os.path.join(source_folder, patient_barcode, 'MuTect2')
		#somaticsniper_folder = os.path.join(source_folder, patient_barcode, 'SomaticSniper')
		#strelka_folder = os.path.join(source_folder, patient_barcode, 'Strelka', 'results')
		#varscan_folder = os.path.join(source_folder, patient_barcode, 'Varscan')

		results = {
			'muse': search_folder(muse_folder, '.Muse.vcf'),
			'mutect': search_folder(mutect_folder, '.mutect2.vcf'),
			'somaticsniper': search_folder(somaticsniper_folder, ".somaticsniper.hq.vcf"),
			'strelka_indel': search_folder(strelka_folder, ".passed.somatic.indels.vcf.strelka.vcf"),
			'strelka-snv': search_folder(strelka_folder, ".passed.somatic.snvs.vcf.strelka.vcf"),
			'varscan-indel': search_folder(varscan_folder, ".raw.indel.vcf"),
			'varscan-snv': search_folder(varscan_folder, ".raw.snp.Somatic.hc.vcf")
		}

		muse_folder, mutect_folder, somaticsniper_folder, strelka_folder, varscan_folder = get_folders(destination_folder, patient_barcode)
		output = list()
		for filenames in results.values():
			if filenames is None: continue
			for source in filenames:
				caller = os.path.split(os.path.split(source)[0])[1]
				basename = os.path.basename(source)
				destination = os.path.join(destination_folder, patient_barcode, caller, basename)
				#Hack for strelka directory structure
				destination = destination.replace('results', 'Strelka')

				if DEBUG: print(source, "->", destination)
				if not os.path.isdir(os.path.dirname(destination)): os.makedirs(os.path.dirname(destination))
				shutil.copy(source, destination)

if __name__ == "__main__":
	source_folder = "G:\\Genomic_Analysis\\1_input_vcfs\\"
	destination_folder = "G:\\test_folder\\"

	transferVariants(source_folder, destination_folder)






