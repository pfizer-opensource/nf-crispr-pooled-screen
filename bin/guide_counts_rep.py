#!/usr/bin/env python

import argparse
import logging
import pandas as pd
import os
import sys
import yaml

from typing import Dict,List

from pooled_utils import is_csv_file,load_csv,parse_sample_name

# NOTE: This multiple representation is not general.
#           Biological replicates -> rep_[1-X]
#           Cell pellet replicates -> [A-Z]
#           PCR replications -> [1-X] this is the representation notion below. How much data is needed to be statistical significant


logger = logging.getLogger(os.path.basename(__file__))

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate the correspondence between the FastQ Files and Samples' Labels.",
        epilog="Example: python generate_sample_label.py <fastq files directory> <list of sample labels>",
    )
    parser.add_argument(
        "--count_file",
        metavar="COUNT_FILE",
        type=str,
        required=False,
        help="Path to the count file.",
    )
    parser.add_argument(
        "--sample_meta_yaml",
        metavar="SAMPLE_META_YAML",
        type=str,
        required=True,
        help="Yaml File with Sample Metadata.",
    )
    parser.add_argument(
        "--prefix",
        metavar="PREFIX",
        type=str,
        required=False,
        default="",
        help="Filename Prefix.",
    )
    parser.add_argument(
        "--output_file_map",
        metavar="OUTPUT_FILE_MAP",
        type=str,
        required=False,
        default="output_file_map.yml",
        help="Filename Prefix.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)

class Count:
    """
    Class to analayze crispr sgRNA counts file .

    """
    # These are the object attributes
    count_file:str=""
    table_count:pd.DataFrame=None # Pandas DataFrame
    meta_data:Dict[str,str]=None # Meta data from the yaml file
    sub_samples:List[str]=[] # Subsamples found (e.g., *_A[1-2], *_B[1-2], *_C[1-2])

    def __init__(self,count_file:str ="", sample_meta:str = "",delete_column:List[str]=["Undetermined"]) -> None:
        """
            Initialize the count and preform some default cleanse .

            Args:
                count_file (str): File containing count values (returned by MaGeck )

                sample_meta (str): Yaml file containing sample metadata

                delete_column(List[str]): List of columns to delete from the count file table

        """

        """
        Mageck returns a tab separated file, but I am using the csv as parser
        changing the delimiter character to \t (tab)
        """

        delimiter = '\t'
        if is_csv_file(count_file,delimiter):
            self.table_count=load_csv(count_file,delimiter)
            self.remove_column(delete_column)

            with open(sample_meta, 'r') as stream:
                try:
                    sample_meta_info = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(f"Failed to parse yaml file {load_csv} with error",exc)
                    raise
            if 'samples' in sample_meta_info.keys():
                self.meta_data=sample_meta_info['samples']
            else:
                self.meta_data=sample_meta_info

        else:
            logger.error(f"The {count_file} is not a CSV file.")

    def _load_csv(self,csv_file:str,delimiter:str) -> None:
        """
        Generate panda dataframe from csv count file
        Args:
            csv_file(str): File to be loaded. It can be a CSV, or any tabular file

            delimiter(str): Delimiter character that separates the columns

        """
        self.table_count = pd.read_csv(csv_file, delimiter=delimiter)

    def remove_column(self, column_name: List[str]) -> None:
        """
        Remove column which name starts with column_name

        Args:
            column_name(List[str]): List of column (partial) names to be removed from analysis

        Raise:
            logger.warning: If the none of the column names provided was not found
        """
        cols_to_drop = [col for col in self.table_count.columns if col.startswith(tuple(column_name))]
        if cols_to_drop:
            self.table_count = self.table_count.drop(cols_to_drop, axis=1)
        else:
            logger.warning(f"The columns {column_name} was not found in the ")

    def combine_guide_counts(self) -> List[Dict]:
        """
        Combine the values of the A1 and A2, A[1-2] and B[1-2], A[1-2] and
        B[1-2] and C[1-2].

        The filename nomenclature follows the
        calc_guide.max_(len<how many letters to combine>)X.count.txt

        Returns:
            df_dict(Dict[str, pd.DataFrame]): Dictionary of
                Key:representation (e.g. 1000) and
                Values: combined guide counts pandas DataFrame
        """
        df_dict = {}
        sample_meta_dict = {}

        rep_pcr_values = list(set([self.meta_data[key]['aliquot'] for key in self.meta_data.keys()]))
        rep_pcr_values.sort()

        # get increments of PCR replicates -> ['A', 'AB', 'ABC']
        rep_pcr_values_combo = [rep_pcr_values[:i+1] for i in range(len(rep_pcr_values))]

        for sub_sample_ind in rep_pcr_values_combo:
            # Although the 2 first lines doesn't change, we need them for
            # each pandas frames.
            df = pd.DataFrame()
            df['sgRNA'] = self.table_count['sgRNA']
            df['Gene'] = self.table_count['Gene']

            # check which samples has the PCR replicate A, A or B, A or B or C
            select_samples = [key for key,value in self.meta_data.items() if value['aliquot'] in sub_sample_ind]
            # Assemble the regex to match the columns to ensure we are not matching partial names
            regexp_list = [f'^{sample}$' for sample in select_samples ]
            # Get the columns that match the regex
            selected_columns = self.table_count.filter(regex='|'.join(regexp_list)).columns
            table_pcr_subset = self.table_count.loc[:,selected_columns]

            # select select_samples keys from self.meta_data
            sub_samples = {key: self.meta_data[key] for key in select_samples}

            # # group dictionary based on the values of the conditions

            group_dict = {}
            # group the dictionary based on the conditions for each bio replicate
            for key,value in sub_samples.items():
                group_values = value.get('group_rep')
                if group_values in group_dict.keys():
                    group_dict[group_values].append(key)
                else:
                    group_dict[group_values] = [key]

            representation = -9999
            # create a dictionary with the samples used to compute the columns values
            tmp_dict = {}
            for group_rep,sample_names in group_dict.items():

                # Assemble the regex to match the columns to ensure we are not matching partial names
                regexp_list = [f'^{sample}$' for sample in sample_names ]
                # Get the columns that match the regex
                sample_columns = table_pcr_subset.filter(regex='|'.join(regexp_list)).columns

                """
                Sum the columns of the same name and matching the sub_sample_ind
                sample1_rep1_A1 + sample1_rep1_A2 or

                sample1_rep1_A1 + sample1_rep1_A2 +
                sample1_rep1_B1 + sample1_rep1_B2

                sample1_rep1_A1 + sample1_rep1_A2 +
                sample1_rep1_B1 + sample1_rep1_B2 +
                sample1_rep1_C1 + sample1_rep1_C2
                """
                # sum of the sample_columns
                sum_colum=table_pcr_subset.loc[:, sample_columns].sum(axis=1)
                col_name=group_rep
                df[col_name]=sum_colum
                sample_list = sample_columns.tolist()
                pcr_list = list(set([sub_samples[sample]['aliquot'] for sample in sample_names]))
                pcr_rep_val=sum([sub_samples[sample]['representation'] for sample in sample_names])
                pcr_list.sort()
                tmp_dict[col_name] = {
                    'samples':sample_list,
                    'aliquot': pcr_list,
                    'replicate':sub_samples[sample_names[0]]['replicate'],
                    'conditions':sub_samples[sample_names[0]]['conditions'],
                    'group':self.meta_data[sample_names[0]]['group'],
                    'group_rep': group_rep,
                    'sample_rep_val': pcr_rep_val,
                    }
                if 'reference' in self.meta_data[sample_names[0]].keys():
                    tmp_dict[col_name]['reference'] = self.meta_data[sample_names[0]]['reference']
                else:
                    tmp_dict[col_name]['is_ref'] = self.meta_data[sample_names[0]].get('is_ref', 0)
                representation = max(pcr_rep_val,representation)
            sample_meta_dict[representation] = tmp_dict
            df_dict[representation] = df

            del df

        return df_dict,sample_meta_dict

def save_pandas_to_csv(
        data_frame:pd.DataFrame,
        outfile:str,
        delimiter:str="\t") -> None:
    """
    Save a Pandas DataFrame into a CSV file, that can be comma separated,
    or using some other delimiter such tab (\t)

    Args:
        data_frame(pd.DataFrame): List of column (partial) names to be removed from analysis

        outfile(str): Output file name to save the DataFrame

        delimiter(str): Delimiter to save file. Defaults to "\\t" (tab)

    """
    data_frame.to_csv(outfile,sep=delimiter, index=False)

def export_output_file_map(file_label:Dict[str,str],outfile:str) -> None:
    """
    export the map with label:file map to be used downstream. This is useful
    to get a reference to a file without the need to parse the filename again

    Args:
        file_label:Dict[str,str]: Dictionary containing representation and
            combined count file

        outfile(str): Output file name
    """
    with open(outfile,"w+") as out_file_map:
        yaml.dump(file_label,out_file_map)



def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # We might decide to use a list with the count file instead
    if args.count_file is None or not os.path.isfile(args.count_file):
        raise FileNotFoundError("Please provide an existing count file!")

    if args.sample_meta_yaml is None or not os.path.isfile(args.sample_meta_yaml):
        raise FileNotFoundError("Please provide an existing sample metadata file!")

    count = Count(count_file=args.count_file, sample_meta=args.sample_meta_yaml)

    combined_count_dict,meta_data = count.combine_guide_counts()

    output_map_file = {}
    for representation,pd_frame in combined_count_dict.items():
        file_name = f"{args.prefix}max_{representation}X.count.txt"
        save_pandas_to_csv(
            pd_frame,
            file_name,
            delimiter='\t'
            )
        yml = f"{args.prefix}max_{representation}X.count.yml"
        export_output_file_map(meta_data[representation],yml)
        # Store the map representation:file_name basename to be used to export
        # the output_file needed downstream
        output_map_file[representation]= os.path.splitext(file_name)[0]

    export_output_file_map(output_map_file,args.output_file_map)





if __name__ == "__main__":
    sys.exit(main())
