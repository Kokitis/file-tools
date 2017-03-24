import csv
import hashlib
import os
if os.name == 'nt':
	genome_folder = "C:\\Users\\Deitrickc\\Google Drive\\Genomics\\"
else: genome_folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/"
import sys
sys.path.append(genome_folder)

def getAvailableFolders(computer_name = None):
	folders = list()

	if computer_name == 'WD' or computer_name is None: #Should come first in case a BAM file is also located in 'DELL'
		folders += [
			"/media/upmc/WD_Partition_2/Genome_Files/"
		]
	if computer_name == 'DELL' or computer_name is None:
		folders += [
			"/media/upmc/Seagate_Backup_Plus_Drive/Genome_Files/Verified/"
		]
	if computer_name == 'LMD' or computer_name is None:
		pass

	return folders

def getComputerName():
	filename = "/home/upmc/Documents/computer_name.txt"
	with open(filename, 'r') as file1:
		computer_name = file1.read()

	return computer_name

def generateFileMd5(filename, blocksize=2**20):
	""" Generates the md5sum of a file.
		Parameters
		----------
			filename: string [PATH]
				Path to a file.
			blocksize: int
				The maximum amount of memory to use when generating the md5sum.
		Returns
		-------
			md5sum: string
				The md5sum of the passed file.
	"""
	m = hashlib.md5()
	with open( filename , "rb" ) as f:
		while True:
			buf = f.read(blocksize)
			if not buf: break
			m.update( buf )
	return m.hexdigest()

def verifyFileStatus(filename, expected_md5sum = None):
	""" Verifies a single file.
		Parameters
		----------
			filename: string [PATH]
				path to a file to test. The file may not exist.
			expected_md5sum: string; default None
				The expected md5sum to test against. If None, the md5sum will not be generated.
		Returns
		-------
			response: dict<>
				* 'filename': ""
				* 'status': ""
				* 'file exists'
				* 'file md5sum':
				* 'expected md5sum'
	"""
	file_status = os.path.isfile(filename)
	if file_status and expected_md5sum is not None:
		file_md5sum = generateFileMd5(filename)
		md5sum_status = file_md5sum == expected_md5sum
		status = file_exists and md5sum_status
	else: 
		md5sum_status = None
		status = file_status

	response = {	
		'filename': filename,
		'status': status,
		'file status': file_status,
		'md5 status': md5sum_status
	}
	return response

def verify(row, folder):
	file_info = gdc_api.file_api(row['id'])
	filename = os.path.join(folder, row['id'], row['filename'])

	index_id = file_info['basic_info']['index_id']
	index_name = file_info['basic_info']['index_name']
	index_filename = os.path.join(folder, row['id'], index_name)
	
	main_response = verify_file(filename, file_info['basic_info']['md5sum'])
	index_response= verify_file(index_filename, index_md5sum)

	return main_response, index_response

def writeManifestFile(rows, filename):
	if len(rows) == 0: return None
	basic_headers = ["id", "filename", "md5", "size", "state"]
	other_headers = sorted(set(rows[0].keys()) - set(basic_headers))
	headers = basic_headers + other_headers

	with open('file_verification_status.tsv', 'w', newline = "") as file2:
		writer = csv.DictWriter(file2, delimiter = '\t', fieldnames = headers)
		writer.writeheader()
		writer.writerows(files)

def verifyManifest(input_manifest, output_manifest):
	with open(input_manifest, 'r') as inputmanifest:
		rows = list(csv.DictReader(inputmanifest, delimiter = '\t'))

	available_folders = getAvailableFolders()
	lines = list()

	for row in rows:
		for folder in available_folder:
			abs_filename = os.path.join(folder, row['id'], row['filename'])
			response = verifyFile(abs_filename, row['md5sum'])
			if response['status']: break
		else:
			abs_filename = ""
			#print("({0}/{1})".format(index,len(reader)),response['status'],row['id'], flush = True)

		new_line = {
			'id': row['id'],
			'filename': row['filename'],
			'md5': row['md5'],
			'size': row['size'],
			'state': row['state'],
			'barcode': row['barcode'],
			'category': row['category'],
			'experimental_strategy': row['experimental_strategy'],
			'patient': row['patient'],
			'tissue': row['tissue'],
			'verification': response['status']
		}
		lines.append(new_line)
	lines = sorted(lines, key = lambda s: s['barcode'])

	writeManifestFile(lines, output_manifest)

if __name__ == "__main__":
	manifest_file = ""
