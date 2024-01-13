#!/usr/bin/env python3

"""
Script : dgeapy_analyze.py
Description: Perform Differential Gene Expression Analysis (DGEA) by identifying
the differentially expressed genes from a dataframe. Takes as input a table
in CSV, TSV, or XLSX format containing gene expression data. The script
applies thresholds for adjusted p-values and fold changes to identify
significant gene expression changes. Any  NaN values in the specified columns
are detected and removed from the analysis. Duplicated indexes are
also removed from the analysis by default.

Generates bar plots and volcano plots to visualize the results.

A directory named "dgeapy_output" will be created in the current working
directory including:
    1. The modified dataframe with new columns for fold change
       and gene regulation.
    2. The generated plots. One bar plot and one volcano plot.

Author : Joan Lluis Pons Ramon
Email : joanlluispons@gmail.com
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from dgeapy_utils import (__version__, Color, TermMsg, CustomHelpFormatter, eprint)
from dgeapy_errors import (InvalidArgumentError, FileNotFoundError,
                           MissingArgumentError, UnsupportedFileFormatError,
                           RequiredColumnNotFoundError, NanValuesInIndexColumnError)


# Default values. Can be changed with command line arguments.
defalut_index_column_name = "index"
defalut_log2_fold_change_column_name = "log2_fold_change"
defalut_fold_change_column_name = "fold_change"
defalut_padj_column_name = "padj"
defalut_regulation_column_name = "regulation"

defalut_fold_change_threshold = 1.5
defalut_padj_threshold = 0.05

defalut_plot_formats = ["png"]

defalut_nan_values = ["", "--"]


def setup_parser() -> argparse.ArgumentParser:

    fmt = lambda prog: CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(
            prog="dgeapy.py analyze",
            formatter_class=fmt,
            description="""
Differential Gene Expression Analysis.
Generates tables with differentially expressed plots to visualize the results.
            """,
            usage="dgeapy.py analyze <TABLE> [OPTIONS]",
            epilog="""
Example:
    dgeapy.py analyze data.csv -o results_dir --padj 0.01 --fold-change 1.5
    dgeapy.py analyze data.xlsx -F png --exclude 'GeneX' --nan-values 'NA'

For more information and documentation, visit <https://github.com/jllpons/dgeapy>.
            """,
            )

    parser.add_argument(
            "table",
            metavar="<TABLE>",
            nargs="?",
            type=str,
            help="Path to the gene expression data file (CSV, TSV, XLSX).",
            )
    parser.add_argument(
            "-o", "--output",
            metavar="DIR",
            default=f"{os.getcwd()}/analyze_output",
            type=str,
            help="Specify the output directory (default: cwd).",
            )
    parser.add_argument(
            "-p", "--padj",
            metavar="FLOAT",
            default=defalut_padj_threshold,
            type=float,
            help=f"Adjusted p-value threshold for significance (default: {defalut_padj_threshold}).",
            )
    parser.add_argument(
            "-f", "--fold-change",
            metavar="FLOAT",
            default=defalut_fold_change_threshold,
            type=float,
            help=f"Fold change threshold for significance (default: {defalut_fold_change_threshold}).",
            )
    parser.add_argument(
            "-F", "--formats",
            metavar="STR",
            nargs="?",
            default=defalut_plot_formats,
            type=str,
            action="append",
            help=f"Output formats for plots (e.g. svg, pdf) (default: {defalut_plot_formats}).",
            )
    parser.add_argument(
            "-e", "--exclude",
            metavar="STR",
            nargs="?",
            default=[],
            type=str,
            action="append",
            help="Exclude indexes matching specified patterns.",
            )
    parser.add_argument(
            "-N", "--nan-values",
            metavar="STR",
            nargs="?",
            default=defalut_nan_values,
            action="append",
            type=str,
            help=f"Strings to recognize as NaN (default: {defalut_nan_values}).",
            )
    parser.add_argument(
            "-k", "--keep-duplicated",
            default=False,
            action="store_true",
            help="Keep duplicated indexes (default: False).",
            )
    parser.add_argument(
            "-I", "--index-column",
            metavar="STR",
            default=defalut_index_column_name,
            type=str,
            help=f"Column name for index (default: {defalut_index_column_name}).",
            )
    parser.add_argument(
            "-L", "--log2fc-column",
            metavar="STR",
            default=defalut_log2_fold_change_column_name,
            type=str,
            help=f"Column name for log2 Fold Change (default: {defalut_log2_fold_change_column_name}).",
            )
    parser.add_argument(
            "-P", "--p-column",
            metavar="STR",
            default=defalut_padj_column_name,
            type=str,
            help=f"Column name for adjusted p-values (default: {defalut_padj_column_name}).",
            )

    return parser


def check_args(args: argparse.Namespace) -> None:
    """
    Check if the given command line arguments are valid.

    Parameters:
        args (argparse.Namespace): The command line arguments.

    Raises:
        MissingArgumentError: If the required arguments are not given.
        FileNotFoundError: If the given table does not exist.
        InvalidArgumentError: If the given arguments are not valid.
    """

    if args.table is None:
        raise MissingArgumentError("<TABLE>", closure="\nSee dgeapy.py analyze --help for more information")

    # Check if the given table exists
    table_path = os.path.abspath(args.table)
    if not os.path.isfile(table_path):
        raise FileNotFoundError(file_name=table_path)

    # Check if the given formats are valid
    valid_formats = ["png", "svg", "pdf", "jpg", "jpeg"]
    for format in args.formats:
        if format not in valid_formats:
            raise InvalidArgumentError(invalid_argument=format,
                                       message=f"{TermMsg.ERROR}: Invalid format",
                                       closure=f"\nValid formats are: {Color.YELLOW}{', '.join(valid_formats)}{Color.RESET}")

    # Check if the given fold change threshold is valid
    if args.fold_change < 0:
        raise InvalidArgumentError(invalid_argument=args.fold_change,
                                   message=f"{TermMsg.ERROR}: Fold change threshold must be a positive number",
                                   closure=f"\nSee dgeapy.py analyze --help for more information")

    # Check if the given adjusted p-value threshold is valid
    if args.padj < 0 or args.padj > 1:
        raise InvalidArgumentError(invalid_argument=args.padj,
                                   message=f"{TermMsg.ERROR}: Adjusted p-value threshold must be a number between 0 and 1",
                                   closure=f"\nSee dgeapy.py analyze --help for more information")


def read_and_validate_table(args: argparse.Namespace) -> pd.DataFrame:
    """
    Try to read the input table and perform the required validations:
        - Check if the table format is supported.
        - Check if the table contains the required columns.
        - Check if the index column contains NaN values.
        - Check if the log2 Fold Change and adjusted p-value columns contain NaN values.
        - Remove duplicated indexes if `--keep-duplicated` arg has not been provided.
        - Remove indexes matching the specified patterns.

    Then, add Fold Change and Regulation columns if they are not present in the table,
    as they are required for the analysis.

    Parameters:
        args (argparse.Namespace): The command line arguments.

    Returns:
        pd.DataFrame: The validated table.

    Raises:
        UnsupportedFileFormatError: If the given table format is not supported
                                    or if the table could not be read.
        RequiredColumnNotFoundError: If the table does not contain the required columns.
        NanValuesInIndexColumnError: If the index column contains NaN values.
    """


    def read_table(table_path: str, nan_values: list[str]) -> pd.DataFrame:
        """
        Read the given table file and return a dataframe.

        Parameters:
            table_path (str): The path to the table file.

        Returns:
            pd.DataFrame: The table loaded into a pandas dataframe object.

        Raises:
            UnsupportedFileFormatError: If the given table format is not supported
                                        or if the table could not be read.
        """

        table_path = os.path.abspath(table_path)

        try:
            if table_path.endswith(".csv"):
                table = pd.read_csv(table_path, na_values=nan_values)
            elif table_path.endswith(".tsv"):
                table = pd.read_csv(table_path, sep="\t", na_values=nan_values)
            elif table_path.endswith(".xlsx"):
                table = pd.read_excel(table_path, na_values=nan_values)
            else:
                raise UnsupportedFileFormatError(file_name=table_path)
        except:
            raise UnsupportedFileFormatError(file_name=table_path,
                                             message=f"{TermMsg.ERROR}: Table format is supported, but could not read it",
                                             closure="\nMaybe check if the table extension matches the actual format")
        return table


    def check_column_names(index_column: str,
                           log2fc_column: str,
                           p_column: str,
                           table: pd.DataFrame) -> None:
        """
        Check if the given column names are present in the dataframe.

        Parameters:
            index_column (str): The name of the index column.
            log2fc_column (str): The name of the log2 fold change column.
            p_column (str): The name of the adjusted p-value column.
            table (pd.DataFrame): The dataframe.

        """

        column_names = table.columns.values.tolist()

        if index_column not in column_names:
            raise RequiredColumnNotFoundError(column_name=index_column,
                                              column_type="index",
                                              closure="\nSee dgeapy.py analyze --help for more information")

        if log2fc_column not in column_names:
            raise RequiredColumnNotFoundError(column_name=log2fc_column,
                                              column_type="log2 fold change",
                                              closure="\nSee dgeapy.py analyze --help for more information")

        if p_column not in column_names:
            raise RequiredColumnNotFoundError(column_name=p_column,
                                              column_type="adjusted p-value",
                                              closure="\nSee dgeapy.py analyze --help for more information")


    def find_nan_indexes(table: pd.DataFrame, index_column: str, colum_to_check: str) -> list[str]:
        """
        Find the indices of the rows containing NaN values in the specified columns.

        Parameters:
            table (pd.DataFrame): The dataframe.
            index_column (str): The name of the index column.
            colum_to_check (str): The name of the column to check.

        Returns:
            list[str]: A list of indices containing NaN values.
        """

        nan_indices = []

        nan_rows = table[table[colum_to_check].isna()]
        if not nan_rows.empty:
            nan_indices.extend(table.loc[nan_rows.index, index_column].values.tolist())

        return nan_indices


    def add_fold_change_column(df: pd.DataFrame, log2_fold_change_column_name: str, fold_change_column_name: str) -> pd.DataFrame:
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

        column_names.insert(fold_change_column_position, fold_change_column_name)

        return df.reindex(columns=column_names)


    def add_regulation_column(df: pd.DataFrame, log2_fold_change_column_name: str, regulation_column_name: str) -> pd.DataFrame:
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
        df[regulation_column_name] = "Unaffected"

        df.loc[df[log2_fold_change_column_name] > 0, regulation_column_name] = "Upregulated"
        df.loc[df[log2_fold_change_column_name] < 0, regulation_column_name] = "Downregulated"


        column_names.insert(regulation_column_position, regulation_column_name)

        return df.reindex(columns=column_names)


    # Try to read the table
    try:
        table = read_table(args.table, args.nan_values)
    except UnsupportedFileFormatError as e:
        raise e

    # Check if the table actually contains the required columns
    try:
        check_column_names(args.index_column, args.log2fc_column, args.p_column, table)
    except RequiredColumnNotFoundError as e:
        raise e

    # Check if the index column contains NaN values
    if not table[table[args.index_column].isna()].empty:
        raise NanValuesInIndexColumnError

    # Look for NaN values in log2 Fold Change and p-adj columns
    # Print a warning message if any NaN value is found and remove the rows containing them
    for column in [args.log2fc_column, args.p_column]:
        nan_indexes = find_nan_indexes(table, args.index_column, column)
        if nan_indexes:
            warn_msg = (f"{TermMsg.WARNING}: {Color.YELLOW}{len(nan_indexes)}{Color.RESET} NaN values have "
                        + f"been found in '{Color.YELLOW}{column}{Color.RESET}' column "
                        + "and will be removed from the analysis. "
                        + "Indexes are: "
                        + str(nan_indexes))
            eprint(warn_msg)

            # Remove the indexes that contain NaN values in log2 Fold Change and p-adj columns
            table = table[table[column].notna()]

    # Set the index column as the index of the dataframe
    table = table.set_index(args.index_column)

    # If `--keep-duplicated` arg has not been provided, remove duplicated indexes
    # keeping just one for each
    if args.keep_duplicated is False:
       duplicated_indexes = table[table.index.duplicated(keep=False)]
       if not duplicated_indexes.empty:
           warn_msg = (f"{TermMsg.WARNING}: {Color.YELLOW}{len(duplicated_indexes)}{Color.RESET} duplicated indexes have "
                       + "been found. Only the first copy for each will be kept. "
                       + "Indexes are: "
                       + str(set(duplicated_indexes.index.values.tolist())))
           eprint(warn_msg)

           # Remove duplicated indexes, keeping just one for each
           table = table.loc[~table.index.duplicated(keep="first")]

    if len(args.exclude):
        for pattern in args.exclude:
            has_pattern = table.index.str.contains(pattern)
            info_msg = (f"{TermMsg.INFO}: {has_pattern.sum()} indexes match the pattern" 
                        + f" '{Color.YELLOW}{pattern}{Color.RESET}' and will"
                        + " be removed from the analysis")
            eprint(info_msg)

            # Remove the indexes that match the pattern
            table = table[~table.index.str.contains(pattern)]

    # Add fold change and regulation columns if they are not present in the table
    if defalut_fold_change_column_name not in table.columns.values.tolist():
        table = add_fold_change_column(table, args.log2fc_column, defalut_fold_change_column_name)
    if defalut_regulation_column_name not in table.columns.values.tolist():
        table = add_regulation_column(table, args.log2fc_column, defalut_regulation_column_name)

    return table


def analyze(args: argparse.Namespace, table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Perform the differential gene expression analysis.

    Parameters:
        args (argparse.Namespace): The command line arguments.
        table (pd.DataFrame): The already validated input table.

    Returns:
        tuple: A tuple containing the differentially expressed genes table (pos 0),
               the upregulated genes table (pos 1), and the downregulated genes table (pos 2).
    """

    deg_table = table[(table[args.p_column] < args.padj)
                      & (abs(table[defalut_fold_change_column_name]) >= args.fold_change)]

    upregulated_table = deg_table[deg_table[defalut_regulation_column_name] == "Upregulated"]
    downregulated_table = deg_table[deg_table[defalut_regulation_column_name] == "Downregulated"]

    info_msg = (f"{TermMsg.INFO}: {Color.YELLOW}{len(deg_table)}{Color.RESET} differentially expressed genes found:\n"
                + f"\t- {Color.YELLOW}{len(upregulated_table)}{Color.RESET} upregulated\n"
                + f"\t- {Color.YELLOW}{len(downregulated_table)}{Color.RESET} downregulated")
    eprint(info_msg)

    return deg_table, upregulated_table, downregulated_table


def write_tables_to_output_dir(output_directory: str,
                               input_table_filename: str,
                               deg_table: pd.DataFrame,
                               upregulated_table: pd.DataFrame,
                               downregulated_table: pd.DataFrame) -> None:
    """
    Write the DGE tables in the specified output directory using the same table
    format as the provided input table.

    Parameters:
        output_directory (str): The output directory.
        input_table_filename (str): The name of the input table file.
        deg_table (pd.DataFrame): The differentially expressed genes table.
        upregulated_table (pd.DataFrame): The upregulated genes table.
        downregulated_table (pd.DataFrame): The downregulated genes table.

    Returns:
        None
    """

    match input_table_filename.split(".")[-1]:
        case "csv":
            deg_table.to_csv(f"{output_directory}/deg.csv")
            upregulated_table.to_csv(f"{output_directory}/upregulated.csv")
            downregulated_table.to_csv(f"{output_directory}/downregulated.csv")
        case "tsv":
            deg_table.to_csv(f"{output_directory}/deg.tsv", sep="\t")
            upregulated_table.to_csv(f"{output_directory}/upregulated.tsv", sep="\t")
            downregulated_table.to_csv(f"{output_directory}/downregulated.tsv", sep="\t")
        case "xlsx":
            deg_table.to_excel(f"{output_directory}/deg.xlsx")
            upregulated_table.to_excel(f"{output_directory}/upregulated.xlsx")
            downregulated_table.to_excel(f"{output_directory}/downregulated.xlsx")


def generate_plots(args: argparse.Namespace, table: pd.DataFrame, deg_table: pd.DataFrame) -> None:
    """
    Generate a barplot and a volcano plot to visualize the results.

    Parameters:
        args (args.Namespace): The command line arguments.
        table (pd.DataFrame): The input table.
        deg_table (pd.DataFrame): The differentially expressed genes table.

    Returns:
        None
    """

    def mk_bar_plot(data: pd.DataFrame,
                    output_directory: str,
                    plot_formats: list[str],
                    padj_threshold: str,
                    foldchange_threshold: str) -> None:
        """
        Generate a bar plot with count and percentage annotations.

        Parameters:
            data (pandas.DataFrame): The data for creating the bar plot.
            output_directory (str): The directory to save the generated plot.
            plot_formats (list): A list of plot formats to save the plot in.

        Returns:
            None

        Generates a bar plot with count and percentage annotations.
        """

        plt.figure(figsize=(11,6))

        count_plot = sns.countplot(data=data,
                                   y="significance",
                                   hue="significance",
                                   palette=("silver", "cornflowerblue", "indianred"),
                                   hue_order=["No significant", "Downregulated", "Upregulated"],
                                   order=["Upregulated", "Downregulated", "No significant"],
                                   width=0.5)

        # Adding count and percentage to the barplot
        total = len(data["significance"])  # Total count
        for p in count_plot.patches:
            count = p.get_width()  # Count of each category
            percentage = (count / total) * 100  # Calculate percentage
            count_plot.annotate(f"{count:.0f} ({percentage:.2f}%)", (count, p.get_y() + p.get_height() / 2), ha="left", va="center")

        count_plot.spines["top"].set_visible(False)
        count_plot.spines["right"].set_visible(False)

        plt.xlabel("Differentially expressed genes", size=10)
        plt.ylabel("", size=10)

        title = ("Adjusted p-value < "
                 + str(padj_threshold)
                 + "\nFold Change >= |"
                 + str(foldchange_threshold)
                 + "|")
        count_plot.set_title(title)

        fig = count_plot.get_figure()

        for format in plot_formats:
            plot_name = f"{output_directory}/barplot.{format}"
            fig.savefig(plot_name, format=format, dpi=300)

            if format == "png":
                plot_name = f"{output_directory}/barplot_transparent-bg.{format}"
                fig.savefig(plot_name, format=format, dpi=300, transparent=True)

        plt.close()


    def mk_volcano_plot(data: pd.DataFrame,
                        log2_fold_change_column_name: str,
                        padj_column_name: str,
                        foldchange_threshold: float,
                        padj_threshold: float,
                        output_directory: str,
                        plot_formats: list[str]) -> None:
        """
        Generate a volcano plot with specified thresholds.

        Parameters:
            data (pandas.DataFrame): The data for creating the volcano plot.
            log2_fold_change_column_name (str): The column name for log2 fold change values.
            padj_column_name (str): The column name for adjusted p-values.
            foldchange_threshold (float): The threshold value for fold change.
            padj_threshold (float): The threshold value for adjusted p-values.
            output_directory (str): The directory to save the generated plot.
            plot_formats (list): A list of plot formats to save the plot in.

        Returns:
            None

        Generates a volcano plot with specified thresholds.
        """

        # Calculating the log2 threshold for FoldChange values
        log2FoldChange_threshold = np.log2(foldchange_threshold)
        # Same for -log10 threshold for padj values
        log10_padj_threshold = -np.log10(padj_threshold)

        # Adding -log10(padj) column to the dataframe
        data["-log10(padj)"] = -np.log10(data[padj_column_name])

        no_sig_count = np.in1d(data["significance"], "No significant").sum()
        up_count = np.in1d(data["significance"], "Upregulated").sum()
        down_count = np.in1d(data["significance"], "Downregulated").sum()

        # I just think it looks better with the counts on the legend
        no_sig_newname = f"No significant ({no_sig_count})"
        up_newname = f"Upregulated ({up_count})"
        down_newname = f"Downregulated ({down_count})"


        data["significance"] = data["significance"].map({
            "No significant": no_sig_newname,
            "Upregulated": up_newname,
            "Downregulated": down_newname,
            })

        plt.figure(figsize=(6,7))

        # Creating the plot
        volcano = sns.scatterplot(data=data,
                                  s=15,
                                  linewidth=0.1,
                                  x=log2_fold_change_column_name,
                                  y="-log10(padj)",
                                  hue="significance",
                                  hue_order=[down_newname, no_sig_newname, up_newname],
                                  palette=("cornflowerblue", "silver", "indianred"))

        # Add lines that will better ilustrate the choosen thresholds:
        volcano.axhline(log10_padj_threshold, zorder=0, c="grey", lw=1, ls="-." )
        volcano.axvline(log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")
        volcano.axvline(-log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")


        # The following two lines remove the top and right borders of the plot.
        volcano.spines["top"].set_visible(False)
        volcano.spines["right"].set_visible(False)

        # Creating the legend title.
        legend_title = ("Adjusted p-value < "
                        + str(padj_threshold)
                        + "\nFold Change >= |"
                        + str(foldchange_threshold)
                        + "|")
        volcano.legend(title=legend_title)

        # Label for x axis
        plt.xlabel("$log_{2}$ Fold change", size=12)
        # Label for y axis
        plt.ylabel("-$log_{10}$ Adjusted p-value", size=12)

        # Set x-axis limits symmetrically around 0
        x_lims = volcano.get_xlim()
        max_abs_x = max(abs(x_lims[0]), abs(x_lims[1]))
        volcano.set_xlim(-max_abs_x, max_abs_x)


        fig = volcano.get_figure()

        # Create the same plot in each specificed format.
        for format in plot_formats:
            plot_name = f"{output_directory}/volcano.{format}"
            fig.savefig(plot_name, format=format, dpi=300)

            if format == "png":
                plot_name = f"{output_directory}/volcano_transparent-bg.{format}"
                fig.savefig(plot_name, format=format, dpi=300, transparent=True)

        plt.close()


    # Preparing values for data representation:
    # Add a new column to the dataframe with the significance of each gene
    table["significance"] = np.where(table.index.isin(deg_table.index), "Significant", "No significant")
    # Replace the "significance" column values with the "regulation" column values (Upregulated/Downregulated)
    table.loc[table["significance"] == "Significant", "significance"] = table[defalut_regulation_column_name]

    mk_bar_plot(table,
                args.output,
                args.formats,
                args.padj,
                args.fold_change)
    mk_volcano_plot(data=table,
                    log2_fold_change_column_name=args.log2fc_column,
                    padj_column_name=args.p_column,
                    foldchange_threshold=args.fold_change,
                    padj_threshold=args.padj,
                    output_directory=args.output,
                    plot_formats=args.formats,)


def main():

    args = setup_parser().parse_args()

    # Validate the given arguments
    try:
        check_args(args)
    except (MissingArgumentError, FileNotFoundError, InvalidArgumentError) as e:
        eprint(f"{e}")
        sys.exit(1)

    # Attempt to read the table and perform the necessary validations
    try:
        table = read_and_validate_table(args)
    except (UnsupportedFileFormatError, RequiredColumnNotFoundError,
            NanValuesInIndexColumnError) as e:
        eprint(f"{e}")
        sys.exit(1)

    # Perform the differential gene expression analysis
    deg_table, upregulated_table, downregulated_table = analyze(args, table)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    # Save the results in the specified output directory
    write_tables_to_output_dir(args.output, args.table,
                               deg_table, upregulated_table, downregulated_table)

    generate_plots(args, table, deg_table)

    eprint(f"{TermMsg.INFO}: Analysis done. See you again soon!")

    sys.exit(0)


if __name__ == "__main__":
    main()
