#!/usr/bin/env python

import sys
import os
import re
import gzip


def main():
	#parse command line options
	from argparse import ArgumentParser
	
	parser = ArgumentParser(description='Transform a newcpgreport from EMBOSS into a tsv file',prog='template_script')
	parser.add_argument('-i','--input',type=str,required=True,help='full path to newcpgreport, gzip format accepted accepted')
	#parser.add_argument('-o','--outout_directory',type=str,required=False, default="passes_filtering", help='The name of the column in sequencing file indicating if the reads quality passed the filtering or not. As this name will probably evolve in the futur, this option can take the new name in input and look for it in the header of the file.')
	parser.add_argument('--verbose',type=bool,required=False,help='option to print out warnings during execution, Defaults to no-verbose (silent mode)')
	args = parser.parse_args()
	
	cpgistat = []

	with open(args.input, mode="rt") as f:

		for i, line in enumerate(f):
			lstrip = line.rstrip()
			if lstrip.startswith('ID') :
				chrm = re.search('ID   (.+)  .+ BP.', lstrip).group(1)
				#print(chrm)

			elif lstrip.startswith('FT   CpG island') :
				cpg_island = re.search('FT\s+CpG island\s+(\d+)..(\d+)', lstrip)
				#print(f'{chrm}	{cpg_island.group(1)}	{cpg_island.group(2)}')	


				

			elif lstrip.startswith('FT                    ') :
				stat = re.search('FT\s+\/.+=(.+)', lstrip).group(1)
				cpgistat.append(stat)


			if lstrip.startswith('FT                    /ObsExp') :
				print(f'{chrm}\t{cpg_island.group(1)}\t{cpg_island.group(2)}\t{cpgistat[0]}\t{cpgistat[1]}\t{cpgistat[2]}\t{cpgistat[3]}')
				cpgistat = []
			##elif lstrip.startswith('//') :			



if __name__ == "__main__":
    main()




