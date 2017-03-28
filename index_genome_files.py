import csv
import hashlib
import os
genome_folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/"
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
		folder += [
			"/media/upmc/LMD_Storage/Genome_Files/Verified/",
			"/home/upmc/Downloads/Genome_Files/"
		]

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

def writeManifestFile(rows, filename):
	if len(rows) == 0: return None
	basic_headers = ["id", "filename", "md5", "size", "state"]
	other_headers = sorted(set(rows[0].keys()) - set(basic_headers))
	headers = basic_headers + other_headers

	with open('file_verification_status.tsv', 'w', newline = "") as file2:
		writer = csv.DictWriter(file2, delimiter = '\t', fieldnames = headers)
		writer.writeheader()
		writer.writerows(files)

def verifyFromManifest(input_manifest, output_manifest):
	""" Verifies files listed in a amanifest file.
		Parameters
		----------
			input_manifest: string [PATH]
				The manifest file to read from.
			output_manifest: string [PATH]
				The manifest file to write the output to.
	"""
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

		new_line = {
			'id': row['id'],
			'filename': row['filename'],
			'md5': row['md5'],
			'size': row['size'],
			'state': row['state'],
			'barcode': row['barcode'],
			'category': row['category'],
			'patient': row['patient'],
			'tissue': row['tissue'],
			'status': response['status'],
			'md5status': response['md5 status'],
			'file status': response['file status']
		}
		lines.append(new_line)
	lines = sorted(lines, key = lambda s: s['barcode'])

	writeManifestFile(lines, output_manifest)

def verifyFileIntegrity(folder):
	""" Verifies the status of all files within a folder, including the index.
	"""

	pass

if __name__ == "__main__":
	input_manifest = "full_manifest.tsv"
	output_manifest = "full_manifest.{0}.tsv".format(getComputerName())
	verifyFromManifest(input_manifest, output_manifest)
