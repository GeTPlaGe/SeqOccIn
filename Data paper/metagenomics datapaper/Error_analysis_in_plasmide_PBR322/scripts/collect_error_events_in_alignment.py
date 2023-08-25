#!/usr/bin/env python3

"""
This script compute error rate from a BAM/FILE of reads aligned to a reference sequence.

:Example:
python compute_error_rate.py -h
"""

# Metadata
__author__ = 'Mainguy Jean - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '0.1'
__email__ = 'support.bioinfo.genotoul@inra.fr'
__status__ = 'dev'


from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from collections import namedtuple
import logging
import re

from typing import List, Dict, Generator, Pattern
import bamnostic as bs

import numpy as np
import pandas as pd
from collections import defaultdict




def parse_cs_tag(cs_tag: str, cs_pattern: Pattern[str]) -> List[str]:
    """
    Parse a CS tag from a SAM line and extract the aligned sequence events.

    :param cs_tag: The CS tag from a SAM line.
    :param cs_pattern: The compiled regular expression pattern for parsing CS segments.
    :return: A list of parsed aligned sequence events.
    """
    cs_parsed = []
    for cs_segment in re.findall(cs_pattern, cs_tag):
        event_type = cs_segment[0]
        assert event_type in [':', '+', '*', '-']
        if event_type == "*":
            cs_parsed.append('*')
        elif event_type == ":":
            cs_parsed.append(':' * int(cs_segment[1:]))
        elif event_type in ["+", '-']:
            cs_parsed.append(event_type * len(cs_segment[1:]))
    return cs_parsed


def open_sam_bam(file: str) -> Generator[str, None, None]:
    """
    Open a SAM or BAM file and yield its lines.

    :param file: Path to the SAM or BAM file.
    :yield: Lines from the SAM or BAM file.
    """
    if file.endswith('.bam'):
        with open(file, 'rb') as fl:
            for obj_bamline in bs.AlignmentFile(fl):
                yield str(obj_bamline)
    else:
        with open(file, 'r') as fl:
            for line in fl:
                yield line


def parse_sam_file(sam_file):
    """
    Parse a SAM file and yield parsed SAM lines.

    :param sam_file: Path to the SAM file.
    :yield: Parsed SAM lines.
    """
    positional_fields = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]

    for i, line in enumerate(open_sam_bam(sam_file)):
        if i % 1000 == 0:
            logging.info(f'Parsing line {i}...')

        if i == 0 or len(line.split('\t')) != len(fields):
            # Determine additional fields in the SAM file
            additional_fields = [new_field.split(':')[0] for new_field in line.split('\t')[len(positional_fields):]]

            fields = positional_fields + additional_fields
            logging.debug(f'{i}: Adding additional fields {len(additional_fields)}')
            SamLineParser = namedtuple('SamLineParser', fields)

        cleaned_line = [
            ':'.join(val.split(':')[2:]) if val.startswith(field) else val
            for field, val in zip(fields, line.strip().split('\t'))
        ]

        samline = SamLineParser._make(cleaned_line)

        yield samline



def count_aln_event_by_reads(sam_file: str, min_cov: int = 50, min_id: int = 50) -> List[Dict[str, int]]:
    """
    Counts alignment events by reads in a SAM/BAM file.

    :param sam_file: Path to the SAM file containing alignments.
    :param min_cov: Minimum coverage threshold for alignment filtering (default: 50).
    :param min_id: Minimum identity threshold for alignment filtering (default: 50).
    :return: A list of dictionaries containing the count of alignment events for each read.
    """

    read_aln_infos = []  # List to store alignment event counts for each read
    cs_pattern = re.compile(r'(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)')  # Regular expression pattern for parsing CS tag

    for i, samline in enumerate(parse_sam_file(sam_file)):

        try:
            cs_parsed = parse_cs_tag(samline.cs, cs_pattern)  # Parse the CS tag from the SAM line
        except AttributeError: 
            logging.warning(f"No CS information was found in {samline.QNAME}")
            continue
            
        cs_parsed_string = ''.join(cs_parsed)

        substitution_count = cs_parsed_string.count('*')
        insertion_count = cs_parsed_string.count('+')
        deletion_count = cs_parsed_string.count('-')
        match_count = cs_parsed_string.count(':')

        aln_length = len(cs_parsed_string)
        prct_ident = 100 * match_count / aln_length
        query_cov = 100 * (match_count + insertion_count + substitution_count) / len(samline.SEQ)

        if query_cov < min_cov or prct_ident < min_id:
            logging.debug(f"Amplicon {samline.QNAME} with identity {prct_ident} and coverage {query_cov} does not pass the identity ({min_id}) and/or coverage thresholds ({min_cov})")

        read_aln_info = {}
        read_aln_info['read_id'] = samline.QNAME.split(';size')[0]
        read_aln_info['match'] = match_count
        read_aln_info['substitution'] = substitution_count
        read_aln_info['deletion'] = deletion_count
        read_aln_info['insertion'] = insertion_count
        read_aln_info['read_len'] = len(samline.SEQ)
        read_aln_info['aln_len'] = len(cs_parsed_string)
        read_aln_info['identity'] = prct_ident
        read_aln_info['query_coverage'] = query_cov
        
        read_aln_infos.append(read_aln_info)

    return read_aln_infos


def parse_arguments():
    """
    Parse the script arguments.

    :return: Parsed arguments.
    """
    parser = ArgumentParser(description="Compute error rate from a BAM/SAM file.",
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b', '--bam_file', required=True, help='Path to the BAM/SAM file for which to compute the error rate.')

    parser.add_argument('--min_identity', default=50, type=int, help='Minimum identity percentage with the reference sequence to consider and report a sequence.  Below this threshold, sequences are ignored.')
    parser.add_argument('--min_coverage', default=50, type=int, help='Minimum coverage percentage with the reference sequence to consider a sequence. Below this threshold, sequences are ignored.')

    parser.add_argument("-v", "--verbose", help="Increase output verbosity.",
                        action="store_true")

    parser.add_argument("-d", "--debug", help="Enable debug mode for more detailed output.",
                        action="store_true")

    args = parser.parse_args()
    return args


def main():

    args = parse_arguments()
    if args.debug:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info('Mode debug ON')

    elif args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
        logging.info('Mode verbose ON')

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s")

    
    read_aln_infos = count_aln_event_by_reads(args.bam_file)
    df_aln_event = pd.DataFrame(read_aln_infos)
    

    df_aln_event = df_aln_event.set_index('read_id')
    
    # keep only one alignement per read. taking the one with most matches.
    df_aln_event = df_aln_event.sort_values(by=['match'], ascending=False)
    df_aln_event = df_aln_event[~df_aln_event.index.duplicated(keep='first')] 

    # summing all errors in a column
    df_aln_event['all_error'] = df_aln_event[['substitution', "insertion", 'deletion']].sum(axis=1)

    read_count_file = 'reads_events_count.tsv'
    logging.info(f'writing {read_count_file}')
    df_aln_event.to_csv(read_count_file, sep='\t', index=True)


if __name__ == '__main__':
    main()
