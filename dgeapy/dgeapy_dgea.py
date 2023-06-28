#!/usr/bin/env python3

"""
Description
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


# Default column names
defalut_index_column_name = "index"
defalut_log2_fold_change_column_name = "log2_fold_change"
defalut_fold_change_column_name = "fold_change"
defalut_padj_column_name = "padj"
defalut_regulation_column_name = "regulation"

# Default thresholds
defalut_fold_change_threshold = 1.5
defalut_padj_threshold = 0.05

# Defalut plot formats
defalut_plot_formats = ["png"]


def add_fold_change_column(df, log2_fold_change_column_name, fold_change_column_name):
    """
    Adds a "{fold_change_column_name}" column to the dataframe after the "{log2_fold_change_column_name}" one.
    Fold Change values are the result of an exponentiation with the number 2 as
    the base integer and the corresponding "{log2_fold_change_column_name} value as the
    exponent.

    Args:
        df (pd.DataFrame): The input dataframe to which the fold change column will be added.
        log2_fold_change_column_name (str): The name of the log2 fold change column.
        fold_change_column_name (str): The desired name for the fold change column.

    Returns:
        pd.DataFrame: The modified dataframe with the "{fold_change_column_name}" column added.
    """

    column_names = df.columns.values.tolist()

    # fold_change column will be placed after log2_fold_change column
    fold_change_column_position = column_names.index(log2_fold_change_column_name) + 1

    df[fold_change_column_name] = np.power(2, abs(df[log2_fold_change_column_name]))

    try:
        column_names.insert(fold_change_column_position, fold_change_column_name)
        return df.reindex(columns=column_names)

    # Pandas will rise and error when trying to reindex form a list with duplicate labels
    except ValueError:
        return df


def add_regulation_column(df, log2_fold_change_column_name, regulation_column_name):
    """
    Adds a "{regulation_column_name}" column to the dataframe after each
    {log2_fold_change_column_name} one. A gene is considered as upregulated if the
    log2 Fold Change value  is positive, and downregulated if the value is negative.

    Args:
        df (pd.DataFrame): The input dataframe to which the regulation columns will be added.
        log2_fold_change_column_name (str): The desired name for the fold change column.
        regulation_column_name (str): The desired name for the regulation column.

    Returns:
        pd.DataFrame: The modified dataframe with the "{regulation_column_name}" column added.
    """

    column_names = df.columns.values.tolist()

    regulation_column_position = column_names.index(log2_fold_change_column_name) + 2

    # Set the initial value of "regulation" column as "unaffected" for all genes
    df[regulation_column_name] = "unaffected"

    df.loc[df[log2_fold_change_column_name] > 0, regulation_column_name] = "upregulated"
    df.loc[df[log2_fold_change_column_name] < 0, regulation_column_name] = "downregulated"


    try:
        column_names.insert(regulation_column_position, regulation_column_name)
        return df.reindex(columns=column_names)

    # Pandas will rise and error when trying to reindex form a list with duplicate labels
    except ValueError:
        return df


def main():

    description = """
    Differential Gene Expression Analysis. Determine the differentially expressed 
    genes from a dataframe.
    """

    parser = argparse.ArgumentParser(
                        description=description,
                        usage="dgeapy.py dgea TABLE [options]"
                        )

    parser.add_argument(
            "dataframe",
            metavar="TABLE",
            nargs="?",
            default="",
            type=str,
            help="path to the CSV, TSV or XLSX file",
            )
    parser.add_argument(
            "-o", "-output-directory",
            metavar="STR",
            default=f"{os.getcwd()}/dgeapy_output",
            type=str,
            help="output directory [Default: $CWD/dgeapy_output]",
            )
    parser.add_argument(
            "--padj",
            metavar="FLOAT",
            default=defalut_padj_threshold,
            type=float,
            help=f"adjusted p-value threshold. LESS THAN THRESHOLD WILL BE APPLIED [Default: {defalut_padj_threshold}]",
            )
    parser.add_argument(
            "--fc",
            metavar="FLOAT",
            default=defalut_fold_change_threshold,
            type=float,
            help=f"fold change threshold. ABSOLUTE VALUE EQUAL OR MORE THAN THRESHOLD WILL BE APPLIED [Default: {defalut_fold_change_threshold}]",
            )
    parser.add_argument(
            "--formats",
            metavar="STR",
            nargs="?",
            default=defalut_plot_formats,
            type=str,
            action="append",
            help=f"plot formats [Default: {defalut_plot_formats}]",
            )
    parser.add_argument(
            "--exclude",
            metavar="STR",
            nargs="?",
            default=[],
            type=str,
            action="append",
            help="string pattern to match against indexes. Matched indexes are excluded"
            )
    parser.add_argument(
            "--nan-values",
            metavar="STR",
            nargs="?",
            default=["", "--"],
            action="append",
            type=str,
            help='strings to recognize as NaN values. Transcripts with NaN padj or NaN fold change will be excluded [Default: ["", "--"]]'
            )
    parser.add_argument(
            "--keep-duplicated",
            default=False,
            action="store_true",
            help="if passed, keep duplicate index values [Default: False]"
            )
    parser.add_argument(
            "--index-column",
            metavar="STR",
            default=defalut_index_column_name,
            type=str,
            help=f"name of the column that will be used as index [Default: {defalut_index_column_name}]",
            )
    parser.add_argument(
            "--log2fc-column",
            metavar="STR",
            default=defalut_log2_fold_change_column_name,
            type=str,
            help=f"name of the column containing the log2 Fold Change values [Default: {defalut_log2_fold_change_column_name}]",
            )
    parser.add_argument(
            "--padj-column",
            metavar="STR",
            default=defalut_padj_column_name,
            type=str,
            help=f"name of the column containing the p-adjusted values [Default: {defalut_padj_column_name}]",
            )

    args = parser.parse_args()

    if not args.dataframe:
        parser.print_help()
        sys.exit("\n** ERROR: The path to the dataframe file is required **")

    df_path = os.path.abspath(args.dataframe)
    if not os.path.isfile(df_path):
        raise FileNotFoundError(f"Could not find file: {df_path}")

    fold_change_threshold = args.fc
    padj_threshold = args.padj
    plot_formats = args.formats

    index_column_name = args.index_column
    log2_fold_change_column_name = args.log2fc_column
    padj_column_name = args.padj_column
    fold_change_column_name = defalut_fold_change_column_name
    regulation_column_name = defalut_regulation_column_name

    output_directory = args.o


    try:
        if df_path.endswith(".csv"):
            df = pd.read_csv(
                        df_path,
                        na_values=args.nan_values,
                        )
        elif df_path.endswith(".tsv"):
            df = pd.read_csv(
                        df_path,
                        na_values=args.nan_values,
                        sep="\t",
                        )
        elif df_path.endswith(".xlsx"):
            df = pd.read_excel(
                        df_path,
                        na_values=args.nan_values,
                        )
        else:
            sys.exit("\n** ERROR: unsupported dataframe format, "
                    + " supported formats are CSV, TSV and XLSX **")
    except:
        sys.exit(f"\n** ERROR: could not read {df_path} **")


    df_column_names = df.columns.values.tolist()
    if index_column_name not in df_column_names:
        sys.exit(f"\n** ERROR: {index_column_name} column not found in the dataframe **")
    if log2_fold_change_column_name not in df_column_names:
        sys.exit(f"\n** ERROR: {log2_fold_change_column_name} column not found in the dataframe **")
    if padj_column_name not in df_column_names:
        sys.exit(f"\n** ERROR: {padj_column_name} column not found in the dataframe **")


    if not df[df[index_column_name].isna()].empty:
        sys.exit(f"\n** ERROR: NaN values have been found in {index_column_name} column **")

    df = df.set_index(index_column_name)

    # Check duplicated values if "--keep-duplicated" arg has not been given
    if args.keep_duplicated is False:
       duplicated_indexes = df[df.index.duplicated(keep=False)]

       if not duplicated_indexes.empty:
           duplicated_indexes_listed = set(duplicated_indexes.index.values.tolist())
           print(
            f"\n** INFO: {str(len(duplicated_indexes_listed))} duplicated indexes have been found. The analysis will keep the first copy for each **"
            + "\nIndexes are: "
            + str(duplicated_indexes.index.values.tolist())
            )
           # Remove duplicated indexes, keeping just one for each
           # NOTE: In pandas, logical NOT is represented by "~"
           df = df.loc[~df.index.duplicated(keep="first")]


    # Removing NaN values
    df = df[df[log2_fold_change_column_name].notna()]
    df = df[df[padj_column_name].notna()]

    # If we need to exclude some index/gene names:
    if len(args.exclude):
        patterns_to_exclude = args.exclude
        for pattern in patterns_to_exclude:
            df = df[~df.index.str.contains(pattern)]

    # Trying to not duplicate columns
    if fold_change_column_name not in df_column_names:
        df = add_fold_change_column(df, log2_fold_change_column_name, fold_change_column_name)
    if regulation_column_name not in df_column_names:
        df = add_regulation_column(df, log2_fold_change_column_name, regulation_column_name)

    # This is the actual filtering for the DEG
    diff_expressed_genes_df = df[
            (df[padj_column_name] < padj_threshold)
            & (abs(df[fold_change_column_name]) >= fold_change_threshold)
            ]

    upregulated_genes_df = diff_expressed_genes_df[diff_expressed_genes_df[regulation_column_name] == "upregulated"]
    downregulated_genes_df = diff_expressed_genes_df[diff_expressed_genes_df[regulation_column_name] == "downregulated"]

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    diff_expressed_genes_df.to_csv(f"{output_directory}/deg.tsv", sep="\t")
    upregulated_genes_df.to_csv(f"{output_directory}/upregulated.tsv", sep="\t")
    downregulated_genes_df.to_csv(f"{output_directory}/downregulated.tsv", sep="\t")

    diff_expressed_genes_df.to_excel(f"{output_directory}/DEG.xlsx")
    upregulated_genes_df.to_excel(f"{output_directory}/upregulated.xlsx")
    downregulated_genes_df.to_excel(f"{output_directory}/downregulated.xlsx")


if __name__ == "__main__":
    main()

