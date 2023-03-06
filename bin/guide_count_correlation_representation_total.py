#!/usr/bin/env python

import argparse
import logging
import pandas as pd
import os
import sys

from typing import Dict,List

from pooled_utils import load_csv,sort_alphanumeric

# NOTE: This multiple representation is not general.
#           Biological replicates -> rep_[1-X]
#           Cell pellet replicates -> [A-Z]
#           PCR replications -> [1-X] this is the representation notion below. How much data is needed to be statistical significant


logger = logging.getLogger(os.path.basename(__file__))

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine the Correlation Tables for all Representation values.",
        epilog="Example: python guide_counts_correlation_representation_total.py --corr_file_list \"correlation_1000x.tsv,correlation_1000x.tsv\" --out_file combined_correlation_table.tsv",
    )
    parser.add_argument(
        "--correlation_file_list",
        metavar="CORR_FILE_LIST",
        type=str,
        required=True,
        help="Path to the count file.",
    )
    parser.add_argument(
        "--out_file",
        metavar="OUT_FILE",
        type=str,
        required=True,
        help="Yaml File with Sample Metadata.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)



def main(argv=None):
    """Coordinate argument parsing and program execution."""
    # compile the regular expression for natural sorting to find numbers in strings


    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # We might decide to use a list with the count file instead
    if args.correlation_file_list is None:
        raise FileNotFoundError("Please provide an existing count file!")

    if args.out_file is None:
        raise FileNotFoundError("Please provide an existing sample metadata file!")

    corr_file_list=[corr for corr in args.correlation_file_list.split(",")]

    # create a dictionary where the keys are the basename of the file and the
    # values are the full path . We need to do this since the hash of the
    # folders can impact the sorting of the files
    file_dict = dict([ (os.path.basename(corr_file),corr_file) for corr_file in corr_file_list])
    # Use the keys to sort the filenames
    corr_file_list = sort_alphanumeric(file_dict.keys())
    # Use the sorted keys to get the full path of the files
    corr_file_list = [file_dict[corr_file] for corr_file in corr_file_list]
    # load the tsv files and merge them into a pandas object
    corr_df = pd.DataFrame()
    for corr_file in corr_file_list:
        corr_df = pd.concat([corr_df,load_csv(corr_file,delimiter='\t')])

    # save the pandas object as a tsv file
    corr_df.to_csv(args.out_file,sep='\t',index=False)


if __name__ == "__main__":
    sys.exit(main())
