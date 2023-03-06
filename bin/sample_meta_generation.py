#!/usr/bin/env python

import argparse
import logging
import glob
import os
import pandas as pd
import re
import sys
import yaml


logger = logging.getLogger(os.path.basename(__file__))

REQUIRE_COLUMNS = ["sample", "replicate", "aliquot",
                    "representation","group"]
# Same as above, but his columns need to be tested as prefix
REQUIRE_COLUMNS_PREFIX= ["condition_"]

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate the correspondence between the FastQ Files and Samples' Labels.",
        epilog="Example: python generate_sample_label.py <fastq files directory> <list of sample labels>",
    )
    parser.add_argument(
        "--fastq_dir",
        metavar="FASTQ_DIR",
        type=str,
        required=True,
        help="Path to the folder containing fastq files.",
    )
    parser.add_argument(
        "--sample_input_meta_table",
        metavar="SAMPLE_META_TABLE.csv",
        type=str,
        required=True,
        help="Input table with sample meta data.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Yaml output file with label and fastq files.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def get_fastq_list_from_dir(fastq_dir):
    """Find the fastq files within a particular folder"""
    fastq_list = [
        os.path.relpath(x) for x in glob.glob(f'{fastq_dir}/*.fastq*') if os.path.isfile(x)
    ]
    return fastq_list

# import csv file and parse with pandas
# docs string in google format
def build_metadata(my_fastq_list,csv_file):
    """Load the sample meta data table csv file into a pandas dataframe"""

    def _define_group(sample):
            # Get the group name based on the conditions
            cond_str=""
            for key,cond in sample['conditions'].items():
                key = key.capitalize()
                cond_str=f"{cond_str}{key} {cond} "

            # remove the last space
            cond_str = cond_str[:-1]
            group = f"{cond_str}"

            return group

    # test if the number of groups are the same as defined by the combination of conditions
    def _test_groups(sample_meta):
        # check if the sample_meta with a particular group has the same conditions
        groups = []
        num_conditions = len(list(sample_meta.values())[0]['conditions'].keys())

        for sample in sample_meta.values():
            if sample['group'] not in groups:
                groups.append(sample['group'])
        # get the conditions for each group
        condition_groups = []
        for group in groups:
            conditions = []
            for sample in sample_meta.values():
                if sample['group'] == group:
                    conditions = conditions + list((sample['conditions'].values()))
            # get the unique conditions
            condition_unique = list(set(conditions))
            #check if condition_unique has more elements than number of conditions keys
            if len(condition_unique) != num_conditions:
                return False
            condition_groups.append(condition_unique)

        # Assuming that one condition per group. If not, return error
        # unique values of list of lists
        unique_condition_groups = list(set(tuple(x) for x in condition_groups))
        if len(unique_condition_groups) != len(groups):
            return False

        return True

    def _extract_group(sample,row):

        # in case group is empty, generate the group name from the other columns
        if not row['group'] or pd.isna(row['group']):
            group  = _define_group(sample)
            row['group'] = group
        # check if group matches the string contains {blah}
        elif re.search(r"{.*}", row['group']):
            # get the values within the all curly brackets
            to_replace = re.findall(r"{(.*?)}", row['group'])
            for value in to_replace:
                # replace the value with the value in the row
                row['group'] = row['group'].replace(
                    f"{{{value}}}", str(row[value])
                )
            group = row['group']
        else:
            # if the group does not contain curly brackets, just return the value
            # store row['group'] into a single entry list
            group = f"{row['group']}"
        return group

    def _define_group_rep(sample,row):
        return f"{row['group']}_rep.{row['replicate']}"
    meta_data = pd.read_csv(csv_file)

    # extract columns names into a list
    columns = meta_data.columns.tolist()

    # check if all the required columns are present and at least on columns
    # with the prefix "condition_"

    if not all([x in columns for x in REQUIRE_COLUMNS]):
        raise Exception(
            f"Missing columns in the input table. Required columns are: {REQUIRE_COLUMNS}"
        )

    for prefix in REQUIRE_COLUMNS_PREFIX:
        if not any([x.startswith(prefix) for x in columns]):
            raise Exception(
                f"Missing columns in the input table. At least one column should have the prefix '{prefix}'"
            )

    # construct a dictionary with main key as sample name and value as a dictionary with the rest of the columns
    sample_meta = {}
    for index, row in meta_data.iterrows():
        sample = row['sample']

        # if the name is empty, set all the values break the loop
        if pd.isna(sample):
            break

        # check if the sample name is unique
        if sample in sample_meta:
            raise Exception(f"Sample name {sample} is not unique.")

        sample_meta[sample] = {}

        sample_meta[sample]['replicate'] = row['replicate']
        sample_meta[sample]['aliquot'] = row['aliquot']
        sample_meta[sample]['representation'] = row['representation']
        # Store any other column value not in REQUIRE_COLUMNS or REQUIRE_COLUMNS_PREFIX
        # Merge the 2 dictionaries.
        sample_meta[sample] =  sample_meta[sample] | dict([(column,row[column])
                                    for column in columns
                                    if column not in REQUIRE_COLUMNS and
                                    not any([column.startswith(prefix)
                                    for prefix in REQUIRE_COLUMNS_PREFIX])]
                                    )
        # set the condition as a dictionary with the columns that have the
        # prefix 'condition_', where key is the column name without prefix and
        # value is the value of the column
        sample_meta[sample]['conditions'] = {
            x.replace('condition_', ''): row[x] for x in columns if x.startswith('condition_')
        }

        sample_meta[sample]['group'] = _extract_group(sample_meta[sample],row)
        sample_meta[sample]['group_rep'] = _define_group_rep(sample_meta[sample],row)

        fastq_files = [x for x in my_fastq_list if sample == parse_fastq_file_name(os.path.basename(x))['name']]

        if len(fastq_files) == 0 :
            raise Exception(f"No fastq files found for sample {sample}")

        for fastq in fastq_files:
            filename_parse = parse_fastq_file_name(os.path.basename(fastq))

            read = f"read{filename_parse['read']}"

            # check if the read is already present in the dictionary
            if read not in sample_meta[sample]:
                sample_meta[sample][read] = []

            sample_meta[sample][read].append(fastq)
        sample_meta[sample]['sequence'] = int(filename_parse['sequence'])

    # test groups
    if not _test_groups(sample_meta):
        raise Exception(f"Number of groups does not match the number of conditions.")

    return sample_meta


def parse_fastq_file_name(fastq):
    pattern = r"(?P<name>[\w\-]+)_(S(?P<sequence>\d+))_(L(?P<lane>\d?)+_)?R?(?P<read>\d)(_\d+)?\.fastq"
    p = re.match(pattern,fastq)

    if p:
        return p.groupdict()
    else:
        logger.error(f"Failed to parse '{fastq}'")
        # raise value error
        raise ValueError(f"Failed to parse '{fastq}'")

def get_metadata_from_dir_fastq_list(fastq_list):
    """Extract samples name list from fastq file path list"""
    sample_dict = [
        parse_fastq_file_name(os.path.basename(fastq))
        for fastq in fastq_list
    ]
    fastq_map = {}
    for fastq_file,pars in zip(fastq_list, sample_dict):
        sample_name = pars["name"]
        read_key = f'read{pars["read"]}'

        # initiate the dictionary for current sample
        if sample_name not in fastq_map:
            fastq_map[sample_name] = {}

        # initiate list for the different lanes if exist
        if read_key not in fastq_map[sample_name]:
            fastq_map[sample_name][read_key] = []

        fastq_map[sample_name]["name"] = sample_name
        fastq_map[sample_name][read_key].append(fastq_file)
        fastq_map[sample_name]["sample_sequence"] = int(pars["ssequence"])

    # Sort nested dictionary by sample_sequence
    ordered_dict = dict(sorted(
        fastq_map.items(),
        key = lambda x: x[1]['sample_sequence']
    ))

    return ordered_dict

def save_output_yaml(my_dict_data, output):
    """Export Dictionary to yaml file containing fastq and label."""

    with open(output, 'w+') as yaml_file:
        yaml.dump(my_dict_data, yaml_file, default_flow_style=False, sort_keys=False)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if not os.path.isdir(args.fastq_dir):

        logger.error(f"Folder {args.fastq_dir} does not exist!")
        raise ValueError(f"Folder {args.fastq_dir} does not exist!")


    if not os.path.isfile(args.sample_input_meta_table):

        logger.error(f"Sample input meta table {args.sample_input_meta_table} does not exist!")
        raise ValueError(f"Sample input meta table {args.sample_input_meta_table} does not exist!")


    my_fastq_list = get_fastq_list_from_dir(args.fastq_dir)

    # if my_fastq_list is empty, raise value error and exit
    if not my_fastq_list:
        logger.error(f"No fastq files found in the directory {args.fastq_dir}.")
        raise ValueError(f"No fastq files found in the directory {args.fastq_dir}.")

    my_sample_table = build_metadata(
                        my_fastq_list,
                        args.sample_input_meta_table
                    )

    # export the dictionary as yaml
    save_output_yaml(my_sample_table, args.output)


if __name__ == "__main__":
    sys.exit(main())
