#!/usr/bin/env python3

"""
Given a list of dataframes, compute all possible intersections of the indexes
of the dataframes. Save one dataframe for each intersection and generate venn
diagrams and upset plots for data visualization.
"""


import argparse
import os
import sys
from typing import Any

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_unweighted, venn3, venn3_unweighted
from upsetplot import from_contents, UpSet
import pandas as pd

from dgeapy_analyze import read_table


__author__ = "Joan Lluis Pons Ramon"
__email__ = "joanlluis@gmail.com"
__license__ = "MIT"


# Defalut values. Can be changed with command line arguments.
defalut_index_column_name = "index"

defalut_plot_formats = ["png"]


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
    Templates for terminal messages.
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
    """
    Setup the argument parser.

    Returns:
    argparse.ArgumentParser: The argument parser.
    """

    fmt = lambda prog: CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(
            prog="dgeapy.py analyze",
            formatter_class=fmt,
            description="""
Computes intersections between multiple data files and generates comprehensive intersection
tables and visualizations.
    """,
    usage="dgeapy.py intersections -f <file1> <file2> [...] -n <name1> <name2> [...] [OPTIONS]",
    epilog="""
examples:
    dgeapy.py intersections -f file1.csv file2.csv -n Experiment1 Experiment2 -o results_dir
    dgeapy.py intersections -f file1.csv file2.csv file3.csv -n Exp1 Exp2 Exp3 --formats png

For more information and documentation, visit <https://github.com/jllpons/dgeapy>
    """,)

    required = parser.add_argument_group("required arguments")
    required.add_argument(
            "-f", "--files",
            metavar="<FILE>",
            nargs="?",
            default=[],
            action="append",
            type=str,
            help="Paths to the data files for intersection analysis.",
            )
    required.add_argument(
            "-n", "--names",
            metavar="STR",
            nargs="?",
            default=[],
            action="append",
            type=str,
            help="Names for the data files to label plots and tables.",
            )

    parser.add_argument(
            "-o", "--output",
            metavar="DIR",
            default=f"{os.getcwd()}/dgeapy_intersections_output",
            type=str,
            help="Specify the output directory for results (default: cwd).",
            )
    parser.add_argument(
            "-i", "--index-column",
            metavar="STR",
            default=defalut_index_column_name,
            type=str,
            help=f"Name of the index column in the data files (default: {defalut_index_column_name}).",
            )
    parser.add_argument(
            "-F", "--formats",
            metavar="STR",
            nargs="?",
            default=defalut_plot_formats,
            type=str,
            action="extend",
            help=f"Output formats for the plots (e.g. svg) (default: {defalut_plot_formats}).",
            )
    parser.add_argument(
            "-N", "--nan-values",
            metavar="STR",
            nargs="?",
            default=["", "--", "NA"],
            action="extend",
            type=str,
            help="Strings to recognize as NaN (default: ['', '--', 'NA']).",
            )
    parser.add_argument(
            "-e", "--exclude",
            metavar="STR",
            nargs="?",
            default=[],
            action="extend",
            type=str,
            help="Exclude indexes matching specified patterns.",
            )

    return parser


def check_args(args: argparse.Namespace) -> None | str:
    """
    Check if the arguments are valid.

    Parameters:
    args (argparse.Namespace): The arguments.

    Returns:
    None | str: None if the arguments are valid, otherwise an error message.
    """

    print(args.files)

    if not args.files:
        return f"{TermMsg.ERROR}: No files were provided"
    elif len(args.files) < 2:
        return f"{TermMsg.ERROR}: At least two files are required"

    if not args.names:
        return f"{TermMsg.ERROR}: No names were provided"
    elif len(args.names) != len(args.files):
        return f"{TermMsg.ERROR}: The number of names must match the number of files"

    for f in args.files:
        if not os.path.isfile(f):
            return f"{TermMsg.ERROR}: Could not find {Color.YELLOW}{f}{Color.RESET}"

    for frmt in args.formats:
        if frmt not in ["png", "svg", "pdf", "jpg", "jpeg"]:
            return (f"{TermMsg.ERROR}: Unsupported format {Color.YELLOW}{frmt}{Color.RESET}"
                    + " (supported formats are png, svg, pdf, jpg and jpeg)")

    return None


def obtain_all_possible_combinations(n: int) -> list[str]:
    """
    Generate all possible combinations of 0s and 1s of length n. Excludes the combination of all 0s.

    Parameters:
    n (int): The desired length of the combinations.

    Returns:
    list: A list of all possible combinations of 0s and 1s, each of length n.

    Example:
    >>> obtain_all_possible_combinations(3)
    ['000', '001', '010', '011', '100', '101', '110', '111']
    """

    combinations = []
    for i in range(2**n):
        # Convert the integer to binary and remove the first two characters (0b)
        combination = bin(i)[2:]
        # Add 0s to the left of the binary number to make it n characters long
        combination = "0"*(n-len(combination)) + combination
        combinations.append(combination)

    # Remove the combination of all 0s
    emtpy_combination = "0"*n
    combinations.remove(emtpy_combination)

    return combinations


def obtain_intersctions(combinations: list[str], dfs: list[pd.DataFrame]) -> dict[str, set]:
    """
    Obtain the intersections of the indexes of the dataframes.

    Parameters:
    combinations (list): A list of all possible combinations of 0s and 1s, each of length n.
    dfs (list): A list of dataframes.

    Returns:
    dict:
        keys (str): The combinations of 0s and 1s.
        values (set): The intersection of the indexes of the dataframes that have a 1 in the combination.
    """

    # Obtain the intersections of the indexes of the dataframes
    intersections = {}

    for combination in combinations:
        indexes = []
        # Get the indexes of the dataframes that have a 1 in the combination
        for i in range(len(combination)):
            # If the combination has a 1 in the i position, add the indexes of the dataframe to the list
            if combination[i] == "1":
                indexes.append(set(dfs[i].index.values.tolist()))

        # If the combination has a 0 in the i position, remove the indexes of the dataframe from the list
        for i in range(len(combination)):
            if combination[i] == "0":
                unwanted_indexes = dfs[i].index.values.tolist()
                # indexes variable is a list of sets, so we need to iterate over it
                for j in range(len(indexes)):
                    for indx in unwanted_indexes:
                        if indx in indexes[j]:
                            indexes[j].discard(indx)

        # Get the intersection of the indexes
        intersection = set.intersection(*indexes)
        if len(intersection):
            # Add the intersection to the dictionary
            intersections[combination] = intersection

    return intersections


def generate_intersections_dfs(intersections: dict[str, set], dfs_dict: dict[str, pd.DataFrame], index_column: str) -> dict[str, pd.DataFrame]:
    """
    Generate the dataframes for each intersection.

    Parameters:
    intersections (dict):
        keys (str): The combinations of 0s and 1s.
        values (set): The intersection of the indexes of the dataframes that have a 1 in the combination.
    dfs_dict (dict):
        keys (str): The names of the dataframes.
        values (pandas.DataFrame): The dataframes.
    Returns:
    dict:
        keys (str): The names of the dataframes.
        values (pandas.DataFrame): The dataframes.
    """

    # Generate the dataframes for each intersection
    intersections_dfs = {}

    for combination, intersection in intersections.items():
        # Get the names of the dataframes that have a 1 in the combination
        dfs_combination = [df_name for i, df_name in enumerate(dfs_dict.keys()) if combination[i] == "1"]

        intersection_df_name = []
        intersection_df = pd.DataFrame()

        for indx, df_name in enumerate(dfs_combination):
            df = dfs_dict[df_name]
            intersection_df_name.append(df_name)
            new_columns = {i: f"{i}_{df_name}" for i in df.columns.values.tolist() if i != index_column}

            if indx == 0:
                df = df.rename(columns=new_columns)
                intersection_df = df[df.index.isin(intersection)]
            else:
                df = df.rename(columns=new_columns)
                intersection_df = intersection_df.merge(df, how="inner", on=index_column)

        intersections_dfs["_and_".join(intersection_df_name)] = intersection_df

    return intersections_dfs


def generate_venn2_diagrams(dfs_dict: dict[str, set], output_directory: str, formats: list[str]) -> None:
    """
    Generate the venn2 diagrams, one with the weights and one without.

    Parameters:
    dfs_dict (dict):
        keys (str): The names of the dataframes.
        values (set): A set containing the indexes of the dataframes.
    output_directory (str): The output directory.
    formats (list): The formats of the plots.

    Returns:
    None
    """

    venn2(dfs_dict.values(), set_labels=dfs_dict.keys())

    for frmt in formats:
        plt.savefig(f"{output_directory}/venn2.{frmt}", dpi=300)

        if frmt == "png":
            plt.savefig(f"{output_directory}/venn2_transparent-bg.{frmt}", transparent=True, dpi=300)

    plt.close()

    venn2_unweighted(dfs_dict.values(), set_labels=dfs_dict.keys())

    for frmt in formats:
        plt.savefig(f"{output_directory}/venn2_unweighted.{frmt}", dpi=300)

        if frmt == "png":
            plt.savefig(f"{output_directory}/venn2_unweighted_transparent-bg.{frmt}", transparent=True, dpi=300)

    plt.close()


def generate_venn3_diagrams(dfs_set_dict: dict[str, set], output_directory: str, formats: list[str]) -> None:
    """
    Generate the venn3 diagrams, one with weights and one without.

    Parameters:
    dfs_dict (dict):
        keys (str): The names of the dataframes.
        values (set): A set containing the indexes of the dataframes.
    output_directory (str): The output directory.
    formats (list): The formats of the plots.

    Returns:
    None
    """

    venn3(dfs_set_dict.values(), set_labels=dfs_set_dict.keys())

    for frmt in formats:
        plt.savefig(f"{output_directory}/venn3.{frmt}", dpi=300)

        if frmt == "png":
            plt.savefig(f"{output_directory}/venn3_transparent-bg.{frmt}", transparent=True, dpi=300)

    plt.close()

    venn3_unweighted(dfs_set_dict.values(), set_labels=dfs_set_dict.keys())

    for frmt in formats:
        plt.savefig(f"{output_directory}/venn3_unweighted.{frmt}", dpi=300)

        if frmt == "png":
            plt.savefig(f"{output_directory}/venn3_unweighted_transparent-bg.{frmt}", transparent=True, dpi=300)

    plt.close()


def generate_upset_plot(dfs_list_dict: dict[str, list[Any]], output_directory: str, formats: list[str]) -> None:
    """
    Generate the upset plot.

    Parameters:
    dfs_list_dict (dict):
        keys (str): The names of the dataframes.
        values (list): A list containing the indexes of the dataframes.
    output_directory (str): The output directory.
    formats (list): The formats of the plots.
    """

    data = from_contents(dfs_list_dict)

    plot = UpSet(
            data,
            subset_size="count",
            show_counts="{:d}",
            sort_by="cardinality",
            sort_categories_by="-input",
            include_empty_subsets=True,
            show_percentages=True,
            element_size=55,
            ).plot()

    for frmt in formats:
        plt.savefig(f"{output_directory}/upset.{frmt}", dpi=300)

        if frmt == "png":
            plt.savefig(f"{output_directory}/upset_transparent-bg.{frmt}", transparent=True, dpi=300)

    plt.close()


def main():

    parser = setup_parser()
    args = parser.parse_args()

    error = check_args(args)
    if error:
        print(error, file=sys.stderr)
        sys.exit(1)

    tables = []
    for f in args.files:
        table, error = read_table(f, args.nan_values)
        if error:
            print(error, file=sys.stderr)
            sys.exit(1)
        tables.append(table)

    for idx, tbl in enumerate(tables):
        if not tbl[tbl[args.index_column].isna()].empty:
            print(f"{TermMsg.ERROR}: NaN values have been found in '{Color.YELLOW}{args.index_column}{Color.RESET}' column",
                  file=sys.stderr)
            sys.exit(1)

        if not tbl[args.index_column].is_unique:
            print(f"{TermMsg.ERROR}: Duplicated values have been found in '{Color.YELLOW}{args.index_column}{Color.RESET}' column"
                  + f" from file '{Color.YELLOW}{args.files[idx]}{Color.RESET}'",
                  file=sys.stderr)
            sys.exit(1)

        tbl.set_index(args.index_column, inplace=True)

    if args.exclude:
        for pattern in args.exclude:
            for idx, tbl in enumerate(tables):
                has_pattern = tbl.index.str.contains(pattern)
                print(f"{TermMsg.INFO}: {Color.YELLOW}({has_pattern.sum()}{Color.RESET}) "
                      + f"indexes matching pattern '{Color.YELLOW}{pattern}{Color.RESET}' have been excluded",
                      file=sys.stderr)
                tables[idx] = tbl[~has_pattern]

    n_files = len(tables)

    table_dict = {} # Will contain tables (for generating intersections tables)
    venn_set_dict = {} # Will contain sets of indexes (for venn diagrams)
    upset_list_dict = {} # Will contain lists of indexes (for upset plots)
    for name, tbl in zip(args.names, tables):
        table_dict[name] = tbl
        venn_set_dict[name] = set(tbl.index.values.tolist())
        upset_list_dict[name] = tbl.index.values.tolist()

    intersections = obtain_intersctions(
                        combinations=obtain_all_possible_combinations(n_files),
                        dfs=tables,)

    # Generating intersections dataframes
    intersections_dfs = generate_intersections_dfs(intersections, table_dict, args.index_column)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    for key, value in intersections_dfs.items():
        key = key.replace(" ", "")
        value.to_csv(f"{args.output}/{key}.tsv", sep="\t")
        value.to_excel(f"{args.output}/{key}.xlsx")

    if n_files == 2:
        generate_venn2_diagrams(venn_set_dict,
                                args.output,
                                args.formats)
        generate_upset_plot(upset_list_dict,
                            args.output,
                            args.formats)

    elif n_files == 3:
        generate_venn3_diagrams(venn_set_dict,
                                args.output,
                                args.formats)
        generate_upset_plot(upset_list_dict,
                            args.output,
                            args.formats)

    else:
        generate_upset_plot(upset_list_dict,
                            args.output,
                            args.formats)


if __name__ == "__main__":
    main()



