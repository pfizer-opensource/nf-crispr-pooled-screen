#!/usr/bin/env python

import csv
import pandas as pd
import re
import warnings

from typing import Dict,List

def load_csv(csv_file:str,delimiter:str) -> pd.DataFrame:
    """
    Generate panda dataframe from csv file
    Args:
        csv_file(str): File to be loaded. It can be a CSV, or any tabular file

        delimiter(str): Delimiter character that separates the columns

    """
    return pd.read_csv(csv_file, delimiter=delimiter)


def is_csv_file(file:str,delimiter:str) -> bool:
    """
    Check if the file is a csv file, even if the extension is not a csv

    Args:
        file(str): File to be tested

        delimiter(str): Character separating the columns

    Returns:
        bool: True: If the file is a CSV or if it can be parse as a tabular
            file, False otherwise

    Raises:
        logger.error: If the file is not a CSV or if it failed to open
    """
    try:
        with open(file,newline='') as csvfile:
            csv.Sniffer().sniff(csvfile.readline(), delimiters=delimiter)
    except csv.Error as error:
        print(error)
        return False
    except FileNotFoundError:
        raise FileNotFoundError(f"The {file} was not found.")

    # Make sure we closed the file before returning
    return True

def parse_sample_name(name:str,pattern:str) -> Dict[str, List[str]]:
    """
    Parse the sample column file names and extract key information:
        name: Samples name, everything until A, B, or C sub samples
        abc: Subsample label (A, B, or C)
        read: read number

    Args:
        generead_name(str): Column name stored in the count file

    Raises:
        logger.error: If it was not possible to parse the column name
    """
    p = re.match(pattern,name)
    if p:
        return p.groupdict()
    else:
        warnings.warn(f"It was not possible to parse {name} with pattern {pattern}")

def sort_alphanumeric( string ):
    """ Sorts the given iterable in the way that is expected.

    Required arguments:
    string -- The iterable to be sorted.

    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(string, key = alphanum_key)
