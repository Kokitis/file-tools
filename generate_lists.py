import os

import gdc_api

def generate_file_md5(filename, blocksize=2**20):
    m = hashlib.md5()
    with open( filename , "rb" ) as f:
        while True:
            buf = f.read(blocksize)
            if not buf: break
            m.update( buf )
    return m.hexdigest()


def verify_file(filename, expected_md5sum, check_md5sum = True):
	""" Verifies a single file
	"""
	file_exists = os.path.isfile(filename)
	if check_md5sum and file_exists:
		file_md5sum = generate_file_md5(filename)
	else:           file_md5sum = None

	if check_md5sum:
		file_is_valid = file_exists and (file_md5sum == expected_md5sum)
	else:
		file_is_valid = file_exists

	response = {	
		'filename': filename,
		'status': file_is_valid,
		'file exists': file_exists,
		'file md5sum': file_md5sum,
		'expected md5sum': expected_md5sum
	}
	return response
def _get_file_ids(filename):
	""" Extracts file ids from the supplied file """
	with open(full_list, 'r') as manifest_file:
		reader = list(csv.DictReader(manifest_file, delimiter = '\t'))

	if 'id' in manifest_file[0].keys():
		id_columns = ['id']
	else:
		id_columns = ['SampleUUID', 'NormalUUID']

	file_ids = list()
	for row in reader:
		file_ids += [f[col] for col in id_columns]
	return file_ids

def _parse_sample_list(sample_list):
	standard_list = list()
	for row in sample_list:

		normal = {
			'id': basic_info['file_id'],
			'filename': basic_info['file_name'],
			'md5': "",
			'size': "",
			'state': "submitted",
			'barcode': row['NormalID'],
			'category': row['category'],
			'experimental_strategy': "",
			'patient': row['PatientID'],
			'tissue': row['tissue'],
			'verification': row['response']['status'],
			'file_location': row['abs_filename'],
		}

def verify_files(reader, folders, computer):
	files = list()
	_num_valid_files = 0
	_num_invalid_files = 0

	print("Checking {0} files...".format(len(reader)))
	for index, file_id in enumerate(reader):
		file_info = gdc_api(file_id, 'files')
		basic_info = file_info['basic_info']
		for folder in folders:
			abs_filename = os.path.join(folder, file_id, basic_info['file_name'])
			response = verify_file(abs_filename, row['md5'], check_md5sum)
			if response['status']: break
		else:
			abs_filename = ""
		
		if response['status']:
			_num_valid_files += 1
		else:
			_num_invalid_files += 1
			abs_filename = ""

		print("({0}/{1})".format(index,len(reader)),response['status'],row['id'], abs_filename,flush = True)

		line = {
			'file_id': file_id,
			'case_id': basic_info['case_id'],
			'file_name': basic_info['file_name'],
			'computer': computer,
			'md5': file_info['md5'],
			'file_size': basic_info['file_size'],
			'barcode': basic_info['sample_barcode'],
			'histology': basic_info['histology'],
			'experimental_strategy': basic_info['experimental_strategy'],
			'patient_barcode': basic_info['patient_barcode'],
			'tissue_type': basic_info['tissue_type'],
			'verification': response['status'],
			'file_location': abs_filename,
			'file_exists': response['file exists']
		}
		files.append(line)
	return files

def export_manifest(files, filename):
	manifest = list()
	#id	filename	md5	size	state	barcode	category	computer	file_exists	file_location	patient	sample type	verification

	for row in files:
		line = {
			'id': row['file_id'],
			'case_id': row['case_id'],
			'filename': row['file_name'],
			'md5': row['md5'],
			'size': row['file_size'],
			'state': 'submitted',
			'barcode': row['barcode'],
			'category': row['histology'],
			'computer': computer,
			'file_location': row['file_location'],
			'patient': row['patient_barcode'],
			'tissue_type': row['tissue_type'],
			'verification': row['verification']
		}
		manifest.append(line)
	fieldnames = ["id"	"filename",	"md5",	"size",	"state",	"barcode",	"category",	"computer",
		"file_location",	"patient",	"sample type",	"verification"]
	
	with open(filename, 'r', newline = '') as file1:
		writer = csv.DictWriter(file1, delimiter = '\t', fieldnames = fieldnames)
		writer.writeheader()
		writer.writerows(manifest)

def export_sample_list(files, filename, computer):
	pass

def annotate_full_list(filename, folder, computer, check_md5sum = False, file_format = 'manifest'):
	"""
		full_list: string [Path]
			The filename of a manifest file or sample list with all samples/files to search for.
			The file must have either an 'id' column or both 'SampleUUID' and 'NormalUUID' columns.
		folders: list<string>
			A list of folders to search through. Assumes that all genomes are saved as [folder]/[id]/[filename].
	"""
	file_ids = _get_file_ids(filename)

	verified_files = verify_files(reader, folders)

	if file_format == 'manifest':
		export_manifest(verified_files, "")
	else:
		export_sample_list(verified_files, "")


source_file = ""
folders = ['']
computer = ""
annotate_full_list()





