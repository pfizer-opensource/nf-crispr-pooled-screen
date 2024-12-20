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
from typing import List,Optional
import yaml

from pooled_utils import is_csv_file,load_csv

logger = logging.getLogger(os.path.basename(__file__))

MANDATORY_ANNOTATION_COLUMNS=["ColumnTarget","ColumnCondition","Default"]
def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Perform QC by comparing guide counts between replicates.",
    )
    parser.add_argument(
        "--count_file",
        metavar="COUNT_FILE",
        type=str,
        required=True,
        help="Guide count table file. (e.g. experiment.max_1000X.count.txt)",
    )
    parser.add_argument(
        "--representation",
        metavar="REPRESENTATION",
        type=int,
        required=True,
        help="Maximum cells per guide representation associated with the count file. (e.g., 1000)",
    )
    parser.add_argument(
        "--samples_yml",
        metavar="SAMPLES_YML",
        type=str,
        required=True,
        help="YAML file with sample metadata. (e.g. experiment.max_1000X.count.yml)",
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
        "--norm_method",
        metavar="NORM_METHOD",
        type=str,
        required=False,
        default="total",
        help="Normalization method (control, median, none, total [default])",
    )
    parser.add_argument(
        "--control_guides",
        metavar="CONTROL_GUIDES_FILE",
        type=str,
        required=False,
        default="",
        help="File containing a list of control guides to use for control normalization",
    )
    parser.add_argument(
        "--pseudocount",
        metavar="NUMBER",
        type=int,
        required=False,
        default=1,
        help="Pseudocount value to use to avoid zero errors with logs (default 1)",
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
        help="The desired log level (default INFO).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )
    return parser.parse_args(argv)

class GuideCounts:
    """
    Class to analyze a guide count table
    """
    # Attributes
    raw_counts:pd.DataFrame = None

    def __init__(self, count_file:str, representation:int, meta_yml:str,
                 base_columns:Optional[List[str]]=None, pseudocount:int=1) -> None:
        delimiter = '\t'
        if is_csv_file(count_file, delimiter):
            self.raw_counts = load_csv(count_file, delimiter)
        else:
            logger.error(f"The {count_file} is not a TSV file.")
            return 1
        self.norm_counts = self.raw_counts.copy()
        self.representation = representation
        self.base_columns = list(base_columns) if base_columns is not None else []
        self.pseudocount = pseudocount
        self._is_control_guide = None

        with open(meta_yml, 'r') as stream:
            try:
                self.sample_meta = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                logger.error(f"Failed to parse yaml file {meta_yml} with error", exc)
                raise

    def annotate_table(self, annotate_dict):

        #Using numpy to make the conditional assignment faster
        result = ['unknown']
        target_column = annotate_dict['ColumnTarget']

        if annotate_dict.get('Default') is not None:
            self.raw_counts[target_column]=annotate_dict["Default"]
            result = [annotate_dict["Default"]]
        else:
            logger.error('Missing the default value for data annotation. '
                         'Using "unknown" as the default.')

        for label, values in annotate_dict.items():
            if label in MANDATORY_ANNOTATION_COLUMNS:
                continue
            else:
                condition = self.raw_counts[annotate_dict["ColumnCondition"]].isin(values)
                result = np.where(condition, label, result)

        if target_column not in self.base_columns:
            self.base_columns.append(target_column)

        self.raw_counts[target_column] = result
        self.norm_counts[target_column] = result

    def load_control_guides(self, control_guide_file:str):

        try:
            with open(control_guide_file, 'r') as f:
                control_guides = [line.strip() for line in f]
                logger.info(f'Loaded {len(control_guides)} control guides.')
        except OSError as e:
            logger.error(f'Could not access control guide file {control_guide_file}: {e}')
            return

        if control_guides:
            self._is_control_guide = self.raw_counts['sgRNA'].isin(control_guides)
            num_controls = self._is_control_guide.sum()
            logger.info(f'Found {num_controls} controls in the guide count table.')

    def normalize_data(self,norm_method:str):

        # Calculate per-sample factors to normalize guide counts to total reads in each
        # sample. Matches implementation of total normalization algorithm from MAGeCK.
        def total_normalization_factors(data:pd.DataFrame):
            avg_reads = data.apply(sum).mean()
            return data.apply(lambda x: avg_reads / x.sum())

        # Sample normalization factors for median normalization, implementing the
        # algorithm used in MAGeCK
        def median_normalization_factors(data:pd.DataFrame):
            # Do not use guides with all zero read counts to calculate size factors
            non_zero_guides = data[data.apply(sum, axis=1) > 0]
            # Calculate guide size factors by finding the geometric mean of each guide
            # across samples and then dividing its count in each sample by its geometric mean
            guide_geo_means = non_zero_guides.apply(lambda x: np.exp(np.mean(np.log(x + 1))),
                                                    axis=1)
            guide_size_factors = non_zero_guides.apply(lambda x: x / guide_geo_means,
                                                       axis=0)
            # The sample size factors should be the median of the guide size factors in
            # each sample. However, MAGeCK does not use the true median but instead uses
            # the following approach, which is also used here for consistence with MAGeCK.
            sample_size_factors = guide_size_factors.apply(lambda x: sorted(x)[len(x) // 2],
                                                           axis=0)
            return 1 / sample_size_factors

        # No normalization, a factor of 1.0 for all samples
        def none_normalization_factors(data:pd.DataFrame):
            return pd.Series([1.0] * data.shape[1], index=data.columns)

        #cast base_columns to tuple
        base_columns=tuple(self.base_columns)

        sample_cols = [col for col in self.raw_counts.columns if not col.startswith(base_columns)]
        counts = self.raw_counts[sample_cols]

        if norm_method == 'control':
            if self._is_control_guide is None or sum(self._is_control_guide) == 0:
                logger.error('No control guides found. Cannot perform control normalization.')
                exit(1)
            # Only use the counts from the control guides to calculate normalization
            # factors, using the median normalization approach
            counts = counts.loc[self._is_control_guide]
            logger.info(f'Using {counts.shape[0]} controls to calculate normalization factors.')
            norm_method = 'median'

        if norm_method == 'median':
            norm_factors = median_normalization_factors(counts)
            logger.debug(f'Initial median normalization factors:\n{norm_factors}')
            inf_factors = norm_factors[norm_factors == float('inf')]
            zero_count_fraction = counts.apply(lambda x: sum(x == 0) / x.size)
            if not inf_factors.empty:
                samples = ", ".join(inf_factors.index.tolist())
                logger.warning('These samples have a median count of zero, so median '
                               f'normalization is not possible (divide by zero):\n{samples}')
                logger.warning('Using total normalization instead.')
                norm_factors = total_normalization_factors(counts)
            elif any(zero_count_fraction > 0.45):
                samples = zero_count_fraction[zero_count_fraction > 0.45]
                logger.warning('The fraction of zero-count guides is >0.45 in these samples'
                               f' and median normalization is unstable:\n{samples}')
                logger.warning('Using total normalization instead.')
                norm_factors = total_normalization_factors(counts)
            elif any(zero_count_fraction > 0.3):
                samples = zero_count_fraction[zero_count_fraction > 0.3]
                logger.warning('The fraction of zero-count guides is >0.3 in these samples'
                               f' and median normalization may be unstable:\n{samples}')
                logger.warning('Consider using total normalization instead.')
        elif norm_method == 'total':
            norm_factors = total_normalization_factors(counts)
        elif norm_method == 'none':
            norm_factors = none_normalization_factors(counts)
        else:
            logger.error(f'Normalization method "{norm_method}" not valid.')
            exit(1)           

        logger.info(f'Final normalization factors:\n{norm_factors}')
        if any(norm_factors < 0.2) or any(norm_factors > 5):
            logger.warning('Some samples have unusually small or large normalization factors'
                           ' (<0.2 or >5). Please review the data normalization and QC.')
        self.norm_counts[sample_cols] = self.raw_counts[sample_cols] * norm_factors

    def long_reshape(self):

        def _get_group(sample):
            # Get the value of the Group for the sample based on sample_meta
            return self.sample_meta[sample]['group']

        def _get_bio_replicate(sample):
            # Get the value of the Replicate for the sample based on sample_meta
            return self.sample_meta[sample]['replicate']

        def _get_representation(sample):
            # Get the value of the Representation for the sample based on sample_meta
            return self.sample_meta[sample]['sample_rep_val']

        # reshape the normalized count table to have the sample name next to the value
        melted = pd.melt(
                self.norm_counts,
                id_vars=self.base_columns,
                var_name="Sample"
            )

        sample_raw = melted['Sample']

        melted['Representation'] = sample_raw.apply(_get_representation)
        melted['Group'] = sample_raw.apply(_get_group)
        melted['Replicate'] = sample_raw.apply(_get_bio_replicate)
        replicates = melted['Replicate'].unique()
        # drop column Sample
        melted.drop("Sample",axis=1,inplace=True)
        cols = [col for col in melted.columns if col != 'Replicate' and col != 'value']

        melted = melted.pivot(index=cols, columns='Replicate', values='value')
        melted = melted.reset_index().rename_axis(None, axis=1).rename_axis(None,axis=0)
        melted = melted.sort_values(by='Group', ignore_index=False)

        # rename replicate columns to Norm.Count_Rep.x x being the number of replicates
        for replica in replicates:
            melted.rename(columns={replica:f'Norm.Count_Rep.{replica}'},inplace=True)

        self.long_shape = melted

    def compute_cors(self):
        """
        Calculate correlation coefficients for each pair of replicates
        """

        def _calc_r_squared(data, repl_pair, log=False):
            col_1, col_2 = repl_pair
            repl_data = data[[col_1, col_2]]
            if log:
                repl_data = np.log(repl_data + self.pseudocount)
            return (repl_data[col_1].corr(repl_data[col_2])) ** 2

        groups = self.long_shape.groupby('Group')

        # get the column names that match Norm.Count_Rep
        replicates = self.long_shape.filter(regex='Norm.Count_Rep').columns.tolist()

        # Pairs of replicates such as 1_2 1_3 2_3
        replicate_pairs = set(itertools.combinations(sorted(replicates), 2))

        replicate_cors = []
        for repl_pair in replicate_pairs:
            replicate_cors.append(
                pd.DataFrame({
                    'Group': groups['Group'].unique().str[0],
                    'Comparison': f"{repl_pair[0]} vs {repl_pair[1]}",
                    'Representation': self.representation,
                    'R-Squared': groups.apply(_calc_r_squared, repl_pair),
                    'Log-Log R-Squared': groups.apply(_calc_r_squared, repl_pair, log=True),
                })
            )

        self.replicate_cors = pd.concat(replicate_cors) if replicate_cors else None

    def plot_replicates(self, prefix):
        """
        Generate plots comparing replicates for each condition
        """

        output_file_base = f"{prefix}max_{self.representation}X.replicates_cor"
        output_file = f"{output_file_base}.pdf"

        # Apply the default theme
        sns.set_theme()
        sns.set_style("whitegrid")

        plot_data = self.long_shape.copy()
        plot_data = plot_data.sort_values('Group')

        # get the column names that match Norm.Count_Rep
        replicates = self.long_shape.filter(regex='Norm.Count_Rep').columns.tolist()
        if len(replicates) == 1:
            logger.warning("Only one replicate was found. No plots were generated.")
            return

        # Add a pseudocount to avoid problems with log/log plotting
        plot_data[replicates] = plot_data[replicates] + self.pseudocount

        # Pairs of replicates such as 1_2 1_3 2_3
        replicate_pairs = set(itertools.combinations(sorted(replicates), 2))
        pdf_list = []
        for repl_pair in replicate_pairs:
            # make the plots in a grid of 2 columns
            grid = sns.FacetGrid(
                plot_data,
                col="Group",
                col_wrap=2,
                sharex= True,
                sharey = True
            )

            # Remove the "{col_var} ="  from the plot title
            grid.set_titles(template='{col_name}')
            grid.figure.suptitle("Correlation Group Comparison")
            # Subtitle with Representation
            grid.figure.text(
                0.5,0.945,
                f' Max Representation {self.representation}',
                ha='center',
                fontsize=12
            )

            # Create the plots where 'Norm.Count_Rep.1' is column name for x axis
            # series and 'Norm.Count_Rep.2' in y axis.
            grid.map(
                sns.regplot,
                repl_pair[0],repl_pair[1],
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
                mask = np.logical_and(
                    self.replicate_cors['Group'] == title,
                    self.replicate_cors['Comparison'] == f"{repl_pair[0]} vs {repl_pair[1]}"
                )
                log_corr_sq = self.replicate_cors.loc[mask, 'Log-Log R-Squared'].values[0]

                # Make it looking nicer
                r_sq_label = f'R^2 = {float(log_corr_sq):.3f}'

                ax.annotate(
                    r_sq_label, xy=(0.05,0.9),
                    bbox=dict(
                        boxstyle="round,pad=0.3",
                        fc="white",
                        ec="grey",
                        lw=1
                    ),
                    fontsize=8,
                    xycoords='axes fraction'
                )
            pdf_list.append(f"{output_file_base}_{repl_pair[0]}_{repl_pair[1]}.pdf")
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

    def export_tables(self, prefix):
        """
        Export normalized counts and long shaped table
        """
        self.norm_counts.round(6).to_csv(
            f'{prefix}max_{self.representation}X.count_normalized.txt',
            sep="\t",
            index=False
        )

        self.long_shape.round(6).to_csv(
            f'{prefix}max_{self.representation}X.guide_data_long.tsv',
            sep="\t",
            index=False
        )

        if self.replicate_cors is not None:
            self.replicate_cors.round(6).to_csv(
                f'{prefix}max_{self.representation}X.replicates_cor.tsv',
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

    base_columns=['sgRNA', 'Gene']

    guide_counts = GuideCounts(args.count_file, args.representation, args.samples_yml,
                               base_columns, args.pseudocount)

    if args.annotate_dict:
        # Annotate table
        annotate_list = define_annotation(args.annotate_dict)
        for annotation in annotate_list:
            guide_counts.annotate_table(annotation)
    if args.control_guides:
        guide_counts.load_control_guides(args.control_guides)

    guide_counts.normalize_data(args.norm_method)
    guide_counts.long_reshape()
    guide_counts.compute_cors()

    guide_counts.plot_replicates(args.prefix)
    guide_counts.export_tables(args.prefix)

if __name__ == "__main__":
    sys.exit(main())
