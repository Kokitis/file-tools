import os
from fuzzywuzzy import fuzz

def search(base_folder, search_term, ignorecase = True, fuzzy = False):
	matches = list()
	for fn in os.listdir(base_folder):
		abs_path = os.path.join(base_folder, fn)
		if os.path.isfile(abs_path):
			if _isMatch(fn, search_term, ignorecase, fuzzy): matches.append(abs_path)
		elif os.path.isdir(abs_path):
			matches += search(abs_path, search_term, ignorecase)

	return matches

def _isMatch(string, term, ignorecase = True, fuzzy = False):
	""" Checks if a string matches a search term.
	"""
	if ignorecase: string = string.lower()

	if not fuzzy:
		is_match = term in string
	else:
		if not isinstance(fuzzy, int): fuzzy = 80
		ratio = fuzz.partial_token_set_ratio(string, term)
		is_match = ratio >= fuzzy

	return is_match

if __name__ == "__main__":
	from pprint import pprint
	search_term = 'batch'
	base_folder = 'C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\'
	matches = search(base_folder, search_term, fuzzy = 80)
	for i in matches:
		print(i)