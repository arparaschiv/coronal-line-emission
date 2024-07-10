#######################################################################
# Document name: getwrd.py
# Created by:    Alin Paraschiv arparaschiv@nso.edu
# Original IDL version: getwrd.pro; R. Sterner 
#
# Change log
# 20190827 ARP ported the IDL version to python
#######################################################################

## The script uses regex( findall() ) to extract words from string 
## Will return the n'th word in a text string s 
## Optional it can return the last m-1 words before n

# CALLING EXAMPLE: python3 getwrd.py s n <or> getwrd.py s n m

def getwrd(*argv):
	import re
## check if input argument string s and position n are given
	if len(argv) <2:
		print("no input string given and/or word number given")
	elif len(argv) ==2:
		s = argv[0]
		n = argv[1]
		res = re.findall(r'\w+', s) 
		return res[n-1]
	elif len(argv) ==3:
		s = argv[0]
		n = argv[1]
		m = argv[2]
		res = re.findall(r'\w+', s) 
		return res[n-m-1:n-1]
	else:
		print("wrong number of input arguments")

## this can make the code run from an exterior environment like bash
import sys
if __name__ == "__main__":
	getwrd(str(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]))