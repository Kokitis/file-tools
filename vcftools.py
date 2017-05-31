import os
import vcf
import shutil

def copyVcf(source, destination):
	with open(source, 'r') as input_file:
		reader = vcf.Reader(input_file)
		if 'Varscan' in source:
			reader.formats['DP4'] = reader.formats['DP4']._replace(num=4)
			reader.formats['DP4'] = reader.formats['DP4']._replace(type='Integer')
		with open(destination, 'w') as file1:
			writer = vcf.Writer(file1, reader)
			for record in reader:
				filterOut = '/' in str(record.ALT[0]) or '/' in record.REF
				if not filterOut:
					try:
						writer.write_record(record)
					except ValueError:
						print(record)
	return destination

def splitVcf(filename, output_folder=None, template = None):
	""" Separates a vcf file into indel and snp sections.
		Parameters
		----------
			filename: path
				THe vcf file to split.
			output_folder: path
				The folder to save the output files to. If None,
				the files will be saved to the folder the original file
				resides.
			output_template: string; default None
				If provided, will be used as the basename for the output files.

		Returns
		-------
			output: dict<>
				* 'indel': path
					The indel file
				* 'snp': path
					The snp file
	"""

	folder, basename = os.path.split(filename)
	basename = os.path.basename(basename)
	snp_filename = os.path.join(output_folder, basename + ".snp.vcf")
	indel_filename=os.path.join(output_folder, basename + ".indel.vcf")
	
	with open(filename) as vcf_file:
		reader 			= vcf.Reader(vcf_file)

		snp_writer 		= vcf.Writer(open(snp_filename,   'w'), reader)
		indel_writer 	= vcf.Writer(open(indel_filename, 'w'), reader)
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

def fixRawCallerOutputs(sample, variants, output_folder, somaticseq_folder):
	"""
		Parameters
	"""
	if isinstance(sample, dict):
		patient_id = sample['PatientID']
	else:
		patient_id = sample

	modify_vjsd_script   = os.path.join(somaticseq_folder, "modify_VJSD.py")
	modify_mutect_script = os.path.join(somaticseq_folder, "modify_Mutect.py")
	
	modified_variants = dict()
	for caller, source in variants.items():
		basename = "{}.{}.corrected.vcf".format(patient_id, caller)
		destination = os.path.join(output_folder, basename)
		if caller == 'varscan':
			#modify VJSD . py −method VarScan2 − i n f i l e input . vcf −o u t f i l e output . vcf
			command = """python3 {program} -method Varscan2 -infile {infile} -outfile {outfile}"""
		elif caller == "somaticsniper":
			command = """python3 {program} -method SomaticSniper -infile {infile} -outfile {outfile}"""
		elif caller == "muse":
			command = """python3 {program} -method MuSE -infile {infile} -outfile {outfile}"""
		else:
			command = None

		if command:
			command = command.format(
				program = modify_vjsd_script,
				infile = filename,
				outfile = destination)
		else:
			shutil.copy2(source, destination)
		modified_variants[caller] = destination

	return modified_variants

if __name__ == "__main__":
	folder = ""

	callset = Callset()
