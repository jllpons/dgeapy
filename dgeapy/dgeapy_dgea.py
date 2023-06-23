#!/usr/bin/env python3

"""
Description
"""

import os
import sys
import argparse

import numpy as np
import pandas as pd


# Column names:
INDEX_COLUMN_NAME = 'index'
LOG2_FOLD_CHANGE_COLUMN_NAME = 'log2_fold_change'
FOLD_CHANGE_COLUMN_NAME = 'fold_change'
PADJ_COLUMN_NAME = 'padj'
REGULATION_COLUMN_NAME = 'regulation'

def add_fold_change_column(df):
    """
    Adds a '{FOLD_CHANGE_COLUMN_NAME}' column to the dataframe after the '{LOG2_FOLD_CHANGE_COLUMN_NAME}' one.
    Fold Change values are the result of an exponentiation with the number 2 as
    the base integer and the corresponding '{LOG2_FOLD_CHANGE_COLUMN_NAME} value as the
    exponent.

    Args:
        df (pd.DataFrame): The input dataframe to which the fold change column will be added.

    Returns:
        pd.DataFrame: The modified dataframe with the '{FOLD_CHANGE_COLUMN_NAME}' column added.
    """

    column_names = df.columns.values.tolist()

    # fold_change column will be placed after log2_fold_change column
    fold_change_column_position = column_names.index(LOG2_FOLD_CHANGE_COLUMN_NAME) + 1

    df[FOLD_CHANGE_COLUMN_NAME] = np.power(2, abs(df[LOG2_FOLD_CHANGE_COLUMN_NAME]))

    column_names = column_names.insert(fold_change_column_position, FOLD_CHANGE_COLUMN_NAME)

    return df.reindex(columns=column_names)


def add_regulation_column(df):
    """
    Adds a '{REGULATION_COLUMN_NAME}' column to the dataframe after each
    {FOLD_CHANGE_COLUMN_NAME} one. A gene is considered as upregulated if the
    log2 Fold Change value  is positive, and downregulated if the value is negative.

    Args:
        df (pd.DataFrame): The input dataframe to which the regulation columns will be added.

    Returns:
        pd.DataFrame: The modified dataframe with the '{REGULATION_COLUMN_NAME}' column added.
    """

    column_names = df.columns.values.tolist()

    regulation_column_position = column_names.index(FOLD_CHANGE_COLUMN_NAME) + 1

    # Set the initial value of 'regulation' column as 'unaffected' for all genes
    df[REGULATION_COLUMN_NAME] = 'unaffected'

    df.loc[df[FOLD_CHANGE_COLUMN_NAME] > 0, REGULATION_COLUMN_NAME] = 'upregulated'
    df.loc[df[FOLD_CHANGE_COLUMN_NAME] < 0, REGULATION_COLUMN_NAME] = 'downregulated'

    column_names = column_names.insert(regulation_column_position, REGULATION_COLUMN_NAME)

    return df.reindex(columns=column_names)


def main():

    description = """
    Differential Gene Expression Analysis. Determine the differentially expressed 
    genes from a dataframe.
    """

    parser = argparse.ArgumentParser(
                        description=description,
                        usage='dgeapy.py dgea CSV [options]'
                        )

    parser.add_argument(
            'dataframe',
            metavar='CSV',
            nargs='?',
            default='',
            type=str,
            help='path to the CSV, TSV or XLSX file',
            )
    parser.add_argument(
            '-o', '-output-directory',
            metavar='FLOAT',
            default=f'{os.getcwd()}/dgeapy_output',
            type=str,
            help='output directory. [Default: CWD/dgeapy_output]',
            )
    parser.add_argument(
            '--padj',
            metavar='FLOAT',
            default=0.05,
            type=float,
            help='adjusted p-value threshold. Less than value will be applied. [Default: 0.05]',
            )
    parser.add_argument(
            '--fc',
            metavar='FLOAT',
            default=1.50,
            type=float,
            help='fold change threshold. Equal or more/less than value will be applied. [Default: 1.50]',
            )
    parser.add_argument(
            '--formats',
            metavar='STR STR',
            nargs='?',
            default=['png'],
            type=str,
            action='append',
            help='plot formats. [Default: png]',
            )
    parser.add_argument(
            '-e', '--exclude',
            metavar='STR STR',
            nargs='?',
            default=[],
            type=str,
            help='string pattern to match against indexes. Matched indexes are excluded'
            )
    parser.add_argument(
            '--nan-values',
            metavar='STR STR',
            nargs='?',
            default=[''],
            type=str,
            help='strings to recognize as NaN values. Transcripts with NaN padj or NaN fold change will be excluded. [Default: ""]'
            )
    parser.add_argument(
            '--keep-duplicated',
            default=False,
            type=bool,
            help='if True, keep duplicate index values. [Default: False]'
            )

    args = parser.parse_args()

    if not args.dataframe:
        parser.print_help()
        sys.exit('\n** ERROR: The path to the dataframe file is required **')

    df_path = os.path.abspath(args.dataframe)
    if not os.path.isfile(df_path):
        raise FileNotFoundError(f"Could not find file: {df_path}")

    FOLD_CHANGE_THRESHOLD = args.fc
    PADJ_THRESHOLD = args.padj
    PLOT_FORMATS = args.formats

    if df_path.endswith('.csv'):
        df = pd.read_csv(
                    df_path,
                    na_values=args.nan,
                    )
    elif df_path.endswith('.tsv'):
        df = pd.read_csv(
                    df_path,
                    na_values=args.nan,
                    sep='\t',
                    )
    elif df_path.endswith('.xlsx'):
        df = pd.read_excel(
                    df_path,
                    na_values=args.nan,
                    )
    else:
        sys.exit('\n** ERROR: unsupported dataframe format, '
                + ' supported formats are CSV, TSV and XLSX **')

    df_column_names = df.columns.values.tolist()
    if INDEX_COLUMN_NAME not in df_column_names:
        sys.exit(f'\n** ERROR: {INDEX_COLUMN_NAME} column not found in the dataframe **')
    if LOG2_FOLD_CHANGE_COLUMN_NAME not in df_column_names:
        sys.exit(f'\n** ERROR: {LOG2_FOLD_CHANGE_COLUMN_NAME} column not found in the dataframe **')
    if FOLD_CHANGE_COLUMN_NAME not in df_column_names:
        sys.exit(f'\n** ERROR: {FOLD_CHANGE_COLUMN_NAME} column not found in the dataframe **')
    if PADJ_COLUMN_NAME not in df_column_names:
        sys.exit(f'\n** ERROR: {PADJ_COLUMN_NAME} column not found in the dataframe **')


    # TODO: check if old index stays in the df
    if not df[df[INDEX_COLUMN_NAME].isna()].empty:
        sys.exit(f'\n** ERROR: NaN values have been found in {INDEX_COLUMN_NAME} column **')
    df = df.set_index(INDEX_COLUMN_NAME)

    # Check duplicated values
    if args.keep is False:
       duplicated_indexes = df[df.index.duplicated(keep=False)]

       if not duplicated_indexes.empty:
           duplicated_indexes_listed = duplicated_indexes.index.values.tolist()
           print(
            f'\n** INFO: {len(duplicated_indexes_listed)} duplicated indexes have been found: \n'
            + duplicated_indexes.index.values.tolist()
            + ' **'
            )
           # Remove duplicated indexes just keeping one for each
           df = df.loc[~df.index.duplicated(keep='first')]


    # Removing NaN values
    df = df[df[LOG2_FOLD_CHANGE_COLUMN_NAME].notna()]
    df = df[df[PADJ_COLUMN_NAME].notna()]

    df = add_fold_change_column(df)
    df = add_regulation_column(df)

    # If we need to exclude some index/gene names:
    if len(args.exclude):
        patterns_to_exclude = args.exclude
        for pattern in patterns_to_exclude:
            # In pandas, logical NOT is '~'
            df = df[~df.index.str.contains(pattern)]




if __name__ == "__main__":
    main()

