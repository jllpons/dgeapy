#!/usr/bin/env python3

"""
Perform Differential Gene Expression Analysis (DGEA) by identifying the
differentially expressed genes from a dataframe. Takes as input a table
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

Usage examples:
    dgeapy.py analyze data.csv -o results_dir --padj 0.01 --fold-change 1.5
    dgeapy.py analyze data.xlsx -F png --exclude 'GeneX' --nan-values 'NA'
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


__author__ = "Joan Lluis Pons Ramon"
__email__ = "joanlluis@gmail.com"
__license__ = "MIT"


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


class Color:
    """
    Color class for terminal output.
    """

    GREEN = "\033[32m"
    RED = "\033[31m"
    YELLOW = "\033[33m"

    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"

    RESET = "\033[0m"

class TermMsg:
    """
    Terminal messages.
    """

    ERROR = f"[dgeapy] {Color.RED}Error{Color.RESET}"
    WARNING = f"[dgeapy] {Color.RED}Warning{Color.RESET}"
    INFO = f"[dgeapy] {Color.GREEN}Info{Color.RESET}"


# <https://stackoverflow.com/questions/18275023/dont-show-long-options-twice-in-print-help-from-argparse>
# I have also added `argparse.RawTextHelpFormatter` to the inheritance
# because I want to be able to use "\n" in the description and epilog.
class CustomHelpFormatter(argparse.RawTextHelpFormatter, argparse.HelpFormatter):
    def __init__(self, prog):
        # Initialize with super from RawTextHelpFormatter and also set our custom widths
        super().__init__(prog, max_help_position=40)

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string


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
            default=f"{os.getcwd()}/dgeapy_output",
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


def check_args(args: argparse.Namespace) -> None | str:
    """
    Check if the given command line arguments are valid.

    Parameters:
        args (argparse.Namespace): The command line arguments.

    Returns:
        None | str: None if the arguments are valid,
                    otherwise a string containing the error message.
    """

    if args.table is None:
        return (f"{TermMsg.ERROR}: No table file given. "
                + f"See '{Color.YELLOW}dgeapy.py analyze --help{Color.RESET}' for more information")

    # Check if the given table exists
    table_path = os.path.abspath(args.table)
    if not os.path.isfile(table_path):
        return (f"{TermMsg.ERROR}: Could not find file "
                + f"'{Color.YELLOW}{table_path}{Color.RESET}'")

    # Check if the given formats are valid
    valid_formats = ["png", "svg", "pdf", "jpg", "jpeg"]
    for format in args.formats:
        if format not in valid_formats:
            return (f"{TermMsg.ERROR}: Invalid format '{Color.YELLOW}{format}{Color.RESET}'. "
                    + f"Valid formats are: {Color.YELLOW}{', '.join(valid_formats)}{Color.RESET}")


    # Check if the given fold change threshold is valid
    if args.fold_change < 0:
        return f"{TermMsg.ERROR}: Fold change threshold must be a positive number"

    # Check if the given adjusted p-value threshold is valid
    if args.padj < 0 or args.padj > 1:
        return f"{TermMsg.ERROR}: Adjusted p-value threshold must be a number between 0 and 1"

    return None


def read_table(table_path: str, nan_values: list[str]) -> tuple[pd.DataFrame | None, None | str]:
    """
    Read the given table file and return a dataframe.

    Parameters:
        table_path (str): The path to the table file.

    Returns:
        pd.DataFrame | str: The dataframe if the table was read successfully,
                            otherwise a string containing the error message.
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
            return None, f"{TermMsg.ERROR}: Unsupported table format '{Color.YELLOW}{table_path}{Color.RESET}'"
    except:
        return None, f"{TermMsg.ERROR}: Table format is supported, but could not read it"

    return table, None


def check_column_names(args: argparse.Namespace, table: pd.DataFrame) -> None | str:
    """
    Check if the given column names are present in the dataframe.

    Parameters:
        args (argparse.Namespace): The command line arguments.
        table (pd.DataFrame): The dataframe.


    Returns:
        None | str: None if the column names are valid,
                    otherwise a string containing the error message.
    """

    column_names = table.columns.values.tolist()

    if args.index_column not in column_names:
        return f"{TermMsg.ERROR}: Index column '{Color.YELLOW}{args.index_column}{Color.RESET}' not found in the table."

    if args.log2fc_column not in column_names:
        return f"{TermMsg.ERROR}: Log2 Fold Change column '{Color.YELLOW}{args.log2fc_column}{Color.RESET}' not found in the table."

    if args.p_column not in column_names:
        return f"{TermMsg.ERROR}: Adjusted p-value column '{Color.YELLOW}{args.p_column}{Color.RESET}' not found in the table."

    return None


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

    count_plot = sns.countplot(
                        data=data,
                        y="significance",
                        hue="significance",
                        palette=("silver", "cornflowerblue", "indianred"),
                        hue_order=["No significant", "Downregulated", "Upregulated"],
                        order=["Upregulated", "Downregulated", "No significant"],
                        width=0.5,
                        )

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

    title = (
            "Adjusted p-value < "
            + str(padj_threshold)
            + "\nFold Change >= |"
            + str(foldchange_threshold)
            + "|"
            )
    count_plot.set_title(title)

    fig = count_plot.get_figure()

    for format in plot_formats:
        plot_name = f"{output_directory}/barplot.{format}"
        fig.savefig(plot_name, format=format, dpi=300)

        if format == "png":
            plot_name = f"{output_directory}/barplot_transparent-bg.{format}"
            fig.savefig(plot_name, format=format, dpi=300, transparent=True)

    plt.close()


def mk_volcano_plot(
        data: pd.DataFrame,
        log2_fold_change_column_name: str,
        padj_column_name: str,
        foldchange_threshold: float,
        padj_threshold: float,
        output_directory: str,
        plot_formats: list[str],
        ):
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
    volcano = sns.scatterplot(
            data=data,
            s=15,
            linewidth=0.1,
            x=log2_fold_change_column_name,
            y="-log10(padj)",
            hue="significance",
            hue_order=[down_newname, no_sig_newname, up_newname],
            palette=("cornflowerblue", "silver", "indianred"),
            )

    # Add lines that will better ilustrate the choosen thresholds:
    # zorder: add to the bottom (all other elements will be added in front)
    # c is for color
    # lw is for line with
    # ls is for line simbol
    volcano.axhline(log10_padj_threshold, zorder=0, c="grey", lw=1, ls="-." )
    volcano.axvline(log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")
    volcano.axvline(-log2FoldChange_threshold, zorder=0, c="grey",lw=1,ls="-.")


    # The following two lines remove the top and right borders of the plot.
    volcano.spines["top"].set_visible(False)
    volcano.spines["right"].set_visible(False)

    # Creating the legend title.
    legend_title = (
            "Adjusted p-value < "
            + str(padj_threshold)
            + "\nFold Change >= |"
            + str(foldchange_threshold)
            + "|"
            )
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


def main():

    args = setup_parser().parse_args()

    error = check_args(args)
    if error:
        print(error, file=sys.stderr)
        sys.exit(1)

    table, error  = read_table(args.table, args.nan_values)
    if error:
        print(error, file=sys.stderr)
        sys.exit(1)

    error = check_column_names(args, table)
    if error:
        print(error, file=sys.stderr)
        sys.exit(1)

    if not table[table[args.index_column].isna()].empty:
        print(f"{TermMsg.ERROR}: NaN values have been found in '{Color.YELLOW}{args.index_column}{Color.RESET}' column",
              file=sys.stderr)

    for column in [args.log2fc_column, args.p_column]:
        nan_indexes = find_nan_indexes(table, args.index_column, column)
        if nan_indexes:
            print(f"{TermMsg.WARNING}: {Color.YELLOW}{len(nan_indexes)}{Color.RESET} NaN values have "
                  + f"been found in '{Color.YELLOW}{column}{Color.RESET}' column "
                  + "and will be removed from the analysis. "
                  + "Indexes are: "
                  + str(nan_indexes),
                  file=sys.stderr)
            table = table[table[column].notna()]

    table = table.set_index(args.index_column)

    # Check duplicated values if "--keep-duplicated" arg has not been given
    if args.keep_duplicated is False:
       duplicated_indexes = table[table.index.duplicated(keep=False)]

       if not duplicated_indexes.empty:
           print(f"{TermMsg.WARNING}: Duplicated indexes have been found. "
                 + "Only the first copy for each will be kept. "
                 + "Indexes are: "
                 + str(set(duplicated_indexes.index.values.tolist())),
                 file=sys.stderr)

           # Remove duplicated indexes, keeping just one for each
           # NOTE: In pandas, logical NOT is represented by "~"
           table = table.loc[~table.index.duplicated(keep="first")]

    if len(args.exclude):
        for pattern in args.exclude:
            has_pattern = table.index.str.contains(pattern)
            print(f"{TermMsg.INFO}: {has_pattern.sum()} indexes match the pattern '{Color.YELLOW}{pattern}{Color.RESET}'"
                  + "\nThese indexes will be removed from the analysis.",
                  file=sys.stderr)
            table = table[~table.index.str.contains(pattern)]

    if defalut_fold_change_column_name not in table.columns.values.tolist():
        table = add_fold_change_column(table, args.log2fc_column, defalut_fold_change_column_name)
    if defalut_regulation_column_name not in table.columns.values.tolist():
        table = add_regulation_column(table, args.log2fc_column, defalut_regulation_column_name)

    # This is the actual filtering for the DEG
    deg_table = table[(table[args.p_column] < args.padj)
                      & (abs(table[defalut_fold_change_column_name]) >= args.fold_change)]

    upregulated_table = deg_table[deg_table[defalut_regulation_column_name] == "Upregulated"]
    downregulated_table = deg_table[deg_table[defalut_regulation_column_name] == "Downregulated"]
    print(f"{TermMsg.INFO}: {Color.YELLOW}{len(deg_table)}{Color.RESET} differentially expressed genes found "
          + f"({Color.YELLOW}{len(upregulated_table)}{Color.RESET} upregulated, "
          + f"{Color.YELLOW}{len(downregulated_table)}{Color.RESET} downregulated).",
          file=sys.stderr)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    if args.table.endswith(".csv"):
        deg_table.to_csv(f"{args.output}/deg.csv")
        upregulated_table.to_csv(f"{args.output}/upregulated.csv")
        downregulated_table.to_csv(f"{args.output}/downregulated.csv")
    elif args.table.endswith(".tsv"):
        deg_table.to_csv(f"{args.output}/deg.tsv", sep="\t")
        upregulated_table.to_csv(f"{args.output}/upregulated.tsv", sep="\t")
        downregulated_table.to_csv(f"{args.output}/downregulated.tsv", sep="\t")
    elif args.table.endswith(".xlsx"):
        deg_table.to_excel(f"{args.output}/deg.xlsx")
        upregulated_table.to_excel(f"{args.output}/upregulated.xlsx")
        downregulated_table.to_excel(f"{args.output}/downregulated.xlsx")

    # Preparing values for data representation
    # 1. Add a new column to the dataframe with the significance of each gene
    table["significance"] = np.where(table.index.isin(deg_table.index), "Significant", "No significant")
    # 2. Replace the "significance" column values with the "regulation" column values (Upregulated/Downregulated)
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

    print(f"{TermMsg.INFO}: Results saved in '{args.output}'",
          file=sys.stderr)
    sys.exit(0)


if __name__ == "__main__":
    main()
