#!/usr/bin/env python
import os
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"

import argparse
import itertools
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pypdf import PdfMerger
import seaborn as sns
import sys
from typing import Dict,List
import yaml

from pooled_utils import is_csv_file,load_csv,sort_alphanumeric

logger = logging.getLogger(os.path.basename(__file__))

MANDATORY_ANNOTATION_COLUMNS=["ColumnTarget","ColumnCondition","Default"]
def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate the correspondence between the FastQ Files and Samples' Labels.",
        epilog="Example: python generate_sample_label.py <fastq files directory> <list of sample labels>",
    )
    parser.add_argument(
        "--repr_list",
        metavar="REP_LIST",
        type=str,
        required=True,
        help="List of representation count files labels (e.g., 1000X).",
    )
    parser.add_argument(
        "--rep_count_files",
        metavar="REP_COUNT_FILES",
        type=str,
        required=True,
        help="List of representation count files. (e.g. calc_rep.max_1000X.count.txt)",
    )
    parser.add_argument(
        "--rep_count_files_yml",
        metavar="REP_COUNT_FILES_YML",
        type=str,
        required=True,
        help="YAML file with the list of representation count files. (e.g. calc_rep.max_1000X.count.yml)",
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
        "--annotate_dict",
        metavar="ANNOTATE_DICT",
        type=str,
        required=False,
        help="YAML file with the annotation dictionary. (e.g. annotate_dict.yml)",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)

class CountRepr:
    """
    Class to analyze the combined count files
    """
    # Attributes
    sample_names:List[str]=[]
    table_count:pd.DataFrame=None

    def __init__(self,repr_label:str,count_file:str,meta_yml:str) -> None:
        # Although this class is very similar to Count in
        # calc_guide_counts_by_representation.py, the 2 codes are
        # independent in their execution and life span, not making much sense
        # creating inherent here.
        delimiter = '\t'
        if is_csv_file(count_file,delimiter):
            self.table_count=load_csv(count_file,delimiter)
        else:
            logger.error(f"The {count_file} is not a CSV file.")
            return 1
        self.repr_label = repr_label
        self.table_count = load_csv(count_file, delimiter="\t")

        with open(meta_yml, 'r') as stream:
            try:
                sample_meta_info = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                logger.error(f"Failed to parse yaml file {load_csv} with error",exc)
                raise
        # Information about the provenience of the data of each column
        self.sample_meta = sample_meta_info

    def annotate_table(self,annotate_dict):

        #Using numpy to make the conditional assignment faster
        result = []

        if annotate_dict["Default"]:
            self.table_count[annotate_dict["ColumnTarget"]]=annotate_dict["Default"]
            result=np.where(True,annotate_dict["Default"],'unknown')
        else:
            logger.error("Missing the default value for data annotation")

        for guide_type,value in annotate_dict.items():

            if guide_type in MANDATORY_ANNOTATION_COLUMNS:
                continue
            else:
                condition=self.table_count[annotate_dict["ColumnCondition"]].isin(value)
                result=np.where(condition,guide_type,result)

        self.table_count[annotate_dict["ColumnTarget"]]=result

    def normalize_data(self,base_columns:list):

        #Normalize guide counts and add pseudocount (1 per million reads) to
        # avoid divide by zero errors
        def normalize_col(col:pd.DataFrame):
            return col.div(col.sum()) * 1000000 + 1

        #cast base_columns to tuple
        base_columns=tuple(base_columns)

        cols = [col for col in self.table_count.columns if not col.startswith(base_columns)]

        self.table_count[cols]=(self.table_count[cols]).apply(normalize_col)

    def long_reshape(self, base_columns:list):

        def _set_group(col_value):
            # get all the keys from conditions for all col_values
            return self.sample_meta[col_value]['group']

        def _set_bio_replicate(col_value):
            # assing the value of the Replicate to the sample row based on sample_meta
            return self.sample_meta[col_value]['replicate']

        def _set_representation(col_value):
            # assing the value of the Replicate to the sample row based on sample_meta
            return self.sample_meta[col_value]['sample_rep_val']


        # reshape the table to have the sample name next to the value
        melted = pd.melt(
                self.table_count,
                id_vars= base_columns,
                var_name="Sample"
            )

        sample_raw=melted['Sample']

        melted['Representation']=sample_raw.apply(_set_representation)
        melted['Group']=sample_raw.apply(_set_group)
        melted['Replicate']=sample_raw.apply(_set_bio_replicate)
        num_replicates = melted['Replicate'].unique()
        # drop column Sample
        melted.drop("Sample",axis=1,inplace=True)
        cols = [col for col in melted.columns if col != 'Replicate' and col != 'value']

        melted = melted.pivot(index=cols, columns='Replicate', values='value')
        melted = melted.reset_index().rename_axis(None, axis=1).rename_axis(None,axis=0)
        melted = melted.sort_values(by='Group', ignore_index=False)

        # rename replicate columns to Norm.Count_Rep.x x being the number of replicates
        for replica in num_replicates:
            melted.rename(columns={replica:f'Norm.Count_Replica.{replica}'},inplace=True)
        # The sorting here is different than the original R script due to the
        # pivot command

        self.long_shape = melted

    def compute_cor(self,representation):
        """
        Calculate correlation coefficient for each pair of replicates
        """
        # columns_to_copy[['Day', 'Condition', 'Representation']]=pd.DataFrame()

        # there must be a better way to define the Day, Condition, and
        # Representation, but at this moment I don't know
        groups = self.long_shape.groupby('Group')
        self.replicate_cor = pd.DataFrame()

        self.replicate_cor['Group']=groups['Group'].unique().str[0]

        # get the column names that match Norm.Count_Replica
        replicas=self.long_shape.filter(regex='Norm.Count_Replica').columns.tolist()
        replicas.sort()

        # combinations of replicates such as combo=1_2 1_3 2_3
        combo_list = list(set(itertools.combinations(sorted(replicas),2)))
        for combo in combo_list:
            combo_label = f"{combo[0]}_{combo[1]}"
            self.replicate_cor[f'Cor Squared {combo_label}'] = groups.apply(lambda x : (x[combo[0]].corr(x[combo[1]])**2))
            self.replicate_cor[f'Logged Cor Squared {combo_label}'] = groups.apply(lambda x : (np.log(x[combo[0]]).corr(np.log(x[combo[1]])))**2)

        self.replicate_cor["Representation"] = representation

    def plot_replicates(self,output_file):
        """
        Generate plots comparing replicates for each condition
        """

        # Apply the default theme
        sns.set_theme()
        sns.set_style("whitegrid")

        temp = self.long_shape
        temp = temp.sort_values('Group')

        # Make a column with the titles of the plots to be used in the grid
        # initiate new_column as dataframe series
        # merge the contents of the columns in the list and separate them with \n
        # temp['Series']=temp[self.cond_types].apply(lambda x: '\n'.join(x.astype(str)), axis=1)


        # get the column names that match Norm.Count_Replica
        replicas=self.long_shape.filter(regex='Norm.Count_Replica').columns.tolist()
        if len(replicas) == 1:
            logger.warning("Only one replicate was found. No plots were generated.")
            return

        # combinations of replicates such as combo=1_2 1_3 2_3
        combo_list = list(set(itertools.combinations(sort_alphanumeric(replicas),2)))
        pdf_list = []
        for combo in combo_list:
            combo_label = f"{combo[0]}_{combo[1]}"
            # make the plots in a grid of 2 columns
            grid = sns.FacetGrid(
                temp,
                col="Group",
                col_wrap=2,
                sharex= True,
                sharey = True
            )

            # Remove the "{col_var} ="  from the plot title
            grid.set_titles(template='{col_name}')
            repr=self.repr_label
            grid.fig.suptitle("Correlation Group Comparison")
            # Subtitle with Representation
            grid.fig.text(
                0.5,0.945,
                f' Max Representation {repr}',
                ha='center',
                fontsize=12
            )


            # Create the plots where 'Norm.Count_Replica.1' is column name for x axis
            # series and 'Norm.Count_Replica.2' in y axis.
            grid.map(
                sns.regplot,
                combo[0],combo[1],
                fit_reg=False,
                color="Black",
                scatter_kws={"s": 5}
            )

            # Make x and y log scale
            plt.xscale("log")
            plt.yscale("log")

            # Add the R^2 = X label, where the value used is the logged correlation
            # squared value computed before
            for ax,title in zip(grid.axes.flat, grid.col_names):
                # Select the Logged Cor Squared value where the conditions match
                self.replicate_cor.to_csv("replicate_cor.csv")
                mask = (self.replicate_cor['Group'] == title)
                log_corr_sq=self.replicate_cor.loc[mask, f'Logged Cor Squared {combo_label}'].values[0]

                # Make it looking nicer
                string = f'R^2 = {"{:.3f}".format(float(log_corr_sq))}'
                ax.annotate(
                    string, xy=(0.05,0.9),
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        fc="white",
                        ec="grey",
                        lw=1
                    ),
                    fontsize=8,
                    xycoords='axes fraction'
                )
            pdf_list.append(output_file.split(".")[0]+"_"+combo_label+".pdf")
            grid.savefig(pdf_list[-1])

        if len(pdf_list) > 1:
            # Save all grids in grid_list to output_file
            # merge the pdfs into one file using pypdf
            merger = PdfMerger()
            for pdf in pdf_list:
                merger.append(pdf)
            merger.write(output_file)
            merger.close()
            # Remove the individual pdfs
            for pdf in pdf_list:
                os.remove(pdf)
        else:
            os.rename(pdf_list[0],output_file)

    def export_tables(self,prefix,repr):
        """
        Export normalized counts and and long shaped table
        """
        self.table_count.round(6).to_csv(
            f'{prefix}normalized_count_representation_{repr}.tsv',
            sep="\t",
            index=False
        )

        self.long_shape.round(6).to_csv(
            f'{prefix}guide_data_long_shape_representation_{repr}.tsv',
            sep="\t",
            index=False
        )

        self.replicate_cor.round(6).to_csv(
            f'{prefix}rep_correlation_representation_{repr}.tsv',
            sep="\t",
            index=False
        )

def define_annotation(annotation_file):
    annotate_dict={}

    # read the yaml file and store in a dictionary
    with open(annotation_file, 'r') as stream:
        try:
            annotate_dict = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            logger.error(f"Failed to parse yaml file {annotation_file} with error",exc)
            raise

    annotate_list = []
    # check if Default, ColumnTarget, ColumnCondition are in the annotation file
    mandatory = MANDATORY_ANNOTATION_COLUMNS
    # Assuming annotation is nested in yaml entry
    for _,annotation in annotate_dict.items():
        for key in mandatory:
            if key not in annotation:
                raise ValueError(f"Annotation file must contain {mandatory} and {key} is missing")
        annotate_list.append(annotation)

    return annotate_list

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    count_rep={}
    repr_list=[rep for rep in (args.repr_list).split(",")]
    count_file=[file for file in (args.rep_count_files).split(",")]
    meta_yml=[file for file in (args.rep_count_files_yml).split(",")]
    base_columns=[str('sgRNA'),str('Gene')]
    prefix=args.prefix

    # Annotate criteria
    annotate_list=[]
    if args.annotate_dict:
        annotate_list=define_annotation(args.annotate_dict)

    # If all the operations fall inside this loop, I would suggest split
    # the representation count file channels in 3 and then run the
    # script per file. Use the channels to parallelize
    for repr,file,yml in zip(repr_list,count_file,meta_yml):
        count_rep[repr] = CountRepr(repr,file,yml)
        if annotate_list:
            # Annotate table (This can be changed via args or read via a file)
            for annotate_dict in annotate_list:
                count_rep[repr].annotate_table(annotate_dict)
                target_column=str(annotate_dict['ColumnTarget'])
                if target_column not in base_columns:
                    base_columns.append(target_column)

        count_rep[repr].table_count.to_csv("table_count.csv")
        # Normalize guide counts to avoid divide by zero errors
        count_rep[repr].normalize_data(base_columns)
        count_rep[repr].long_reshape(base_columns)
        count_rep[repr].compute_cor(repr)
        count_rep[repr].plot_replicates(f"{prefix}replicate_coor_{repr}.pdf")

        count_rep[repr].export_tables(prefix,repr)

if __name__ == "__main__":
    sys.exit(main())




