#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A crude slow all python version of dfits

Usage:

apero_dfits {filename}

apero_dfits '{glob command}' --fitsort {header keys}

i.e.

apero_dfits '*.fits' --fitsort OBJECT "HIERARCH ESO DPR TYPE"


Created on 2023-02-14

@author: artigau, cook
"""
import glob
import sys

from astropy.io import fits
from astropy.table import Table


# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def get_args():
	"""
	Get arguments in a crude way (don't use argparse here)
	:return:
	"""
	arglist = sys.argv[1:]

	# find whether we have the fitsort option
	if '--fitsort' in arglist:
		fitsort_args = arglist[arglist.index('--fitsort') + 1:]
		# push all previous arguments into dfits
		dfits_args = arglist[:arglist.index('--fitsort')]
	else:
		fitsort_args = None
		# all args go into dfits
		dfits_args = arglist
	# return the dfits and fitsort arguments
	return dfits_args, fitsort_args


def get_key(key, hdr):
	if key not in hdr and '.' in key:
		eso_key = 'HIERARCH ESO {0}'.format(key.replace('.', ' '))
		return eso_key
	# otherwise return key
	return key


def main():
	"""
	Main function acts like dfits and fitsort (crude and slow)

	:return:
	"""
	dfits_args, fitsort_args = get_args()
	# get files using glob
	files = glob.glob(' '.join(dfits_args))
	# -------------------------------------------------------------------------
	# deal with printing full header
	if fitsort_args is None:
		# loop around files
		for filename in files:
			hdr = fits.getheader(filename)
			# loop around keys and print
			for key in hdr.keys():
				print(key, hdr[key], f'//{hdr.comments[key]}')
			# blank lines
			print('\n\n')
	# -------------------------------------------------------------------------
	# deal with fitsort arguments
	else:
		# storage dict for columns
		storage = dict(FILENAME=[])
		# loop around files
		for it, filename in enumerate(files):
			# load header
			hdr = fits.getheader(filename)
			# add filename to storage
			storage['FILENAME'].append(filename)
			# loop around fits keys
			for rawkey in fitsort_args:
				# get key to check
				key = get_key(rawkey, hdr)
				# add keys for first extension
				if it == 0:
					storage[key] = []
				# deal with new key
				elif key not in storage:
					for key_jt in storage:
						storage[key_jt] = [''] * (it + 1)
				# add key
				if key not in hdr:
					storage[key].append('')
				else:
					storage[key].append(hdr[key])
		# convert to table
		table = Table(storage)
		# print table
		table.pprint(max_lines=-1, max_width=-1)


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
	# run main with no arguments (get from command line - sys.argv)
	ll = main()

# =============================================================================
# End of code
# =============================================================================
