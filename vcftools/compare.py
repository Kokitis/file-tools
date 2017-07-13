
import vcf
from pprint import pprint

def getVcfPositions(filename, chromosome = None):

	all_positions = list()
	filtered_positions = list()

	with open(filename, 'r') as vcf_file:
		reader = vcf.Reader(vcf_file)
		for record in reader:
			if chromosome is None or record.CHROM == chromosome:
				position = (record.CHROM, record.POS)
				all_positions.append(position)
				if record.FILTER:
					filtered_positions.append(position)

	result = {
		'allPositions': all_positions,
		'filteredPositions': filtered_positions
	}

	return result

def compare(left_vcf, right_vcf, chromosome = None):
	"""
		Parameters
		----------
			left_vcf, right_vcf: string
				Paths t othe vcf files to compare.
			chromosome: str; default None
				If provided, the comparison will only be done over this chromosome.
	"""
	left_positions = getVcfPositions(left_vcf, chromosome)
	right_positions= getVcfPositions(right_vcf, chromosome)

	all_left_positions = set(left_positions['allPositions'])
	all_right_positions= set(right_positions['allPositions'])

	filtered_left_positions = set(left_positions['filteredPositions'])
	filtered_right_positions= set(right_positions['filteredPositions'])

	common_positions= all_left_positions & all_right_positions

	unique_positions_left = sorted(all_left_positions - common_positions)
	unique_positions_right= sorted(all_right_positions - common_positions)

	print("Unique positions in the left file: ")
	print(len(unique_positions_left))
	print(len(filtered_left_positions & set(unique_positions_left)))

	print("Unique positions in the right file: ")
	print(len(unique_positions_right))
	print(len(filtered_right_positions & set(unique_positions_right)))

	total_common_positions = len(common_positions)
	total_unique_positions = len(set(unique_positions_left + unique_positions_right))

	total_positions = len(all_left_positions | all_right_positions)
	ratio = total_common_positions / total_positions

	print("Total Common: ", total_common_positions)
	print("Total Unique: ", total_unique_positions)
	print("Total All: ", total_positions)
	print("Ratio: ", ratio)


if __name__ == "__main__":
	left_vcf = "C:\\Users\\Deitrickc\\Downloads\\Downloads\\8a9e7682-0071-441d-8cc1-3bf4b27305e9\\8a9e7682-0071-441d-8cc1-3bf4b27305e9.vcf"
	right_vcf = "C:\\Users\\Deitrickc\\Documents\\Projects\\TCGA-2H-A9GF-CHR1\\SomaticSniper\\TCGA-2H-A9GF-CHR1-11A_vs_TCGA-2H-A9GF-CHR1-01A.somaticsniper.vcf"

	left_vcf = "C:\\Users\\Deitrickc\\Downloads\\Downloads\\3c1bff3c-2210-4d6f-b898-d9268c6a6924\\3c1bff3c-2210-4d6f-b898-d9268c6a6924.snp.Somatic.hc.vcf"
	right_vcf = "C:\\Users\\Deitrickc\\Documents\\Projects\\TCGA-2H-A9GF-CHR1\\Varscan\\TCGA-2H-A9GF-CHR1-11A_vs_TCGA-2H-A9GF-CHR1-01A.varscan.snp.Somatic.hc.vcf"

	compare(left_vcf, right_vcf, 'chr1')

	
