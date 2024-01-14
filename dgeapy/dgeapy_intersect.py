#!/usr/bin/env python3

"""
Script: dgeapy_intersect.py
Description: Given a list of dataframes, compute all possible intersections of
the indexes of the dataframes. Save one dataframe for each intersection and
generate venn diagrams and upset plots for data visualization.

Author: Joan Lluis Pons Ramon
Email: joanlluispons@gmail.com
"""


import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_unweighted, venn3, venn3_unweighted
import pandas as pd
from upsetplot import from_contents, UpSet

from dgeapy_errors import (MissingArgumentError, InvalidArgumentError,
                           FileNotFoundError, UnsupportedFileFormatError,
                           NanValuesInIndexColumnError, DuplicatesInIndexColumnError)
from dgeapy_utils import (Color, TermMsg, CustomHelpFormatter, eprint)


# Defalut values. Can be changed with command line arguments.
default_index_column_name = "index"

default_plot_formats = ["png"]


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
    usage="dgeapy.py intersections -f <file1> -f <file2> [...] -n <name1> -n <name2> [...] [OPTIONS]",
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
            default=f"{os.getcwd()}/intersect_output",
            type=str,
            help="Specify the output directory for results (default: cwd).",
            )
    parser.add_argument(
            "-i", "--index-column",
            metavar="STR",
            default=default_index_column_name,
            type=str,
            help=f"Name of the index column in the data files (default: {default_index_column_name}).",
            )
    parser.add_argument(
            "-F", "--formats",
            metavar="STR",
            nargs="?",
            default=default_plot_formats,
            type=str,
            action="extend",
            help=f"Output formats for the plots (e.g. svg) (default: {default_plot_formats}).",
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


def check_args(args: argparse.Namespace) -> None:
    """
    Check if the arguments are valid.

    Parameters:
    args (argparse.Namespace): The arguments.

    Returns:
        None

    Raises:
        MissingArgumentError: If a required argument is missing.
        InvalidArgumentError: If an argument is invalid.
        FileNotFoundError: If a file does not exist.
    """

    print(args.files)

    if not args.files:
        raise MissingArgumentError("<TABLE>", closure="\nSee dgeapy.py intersect --help for more information")
    elif len(args.files) < 2:
        raise InvalidArgumentError(invalid_argument="<TABLE>",
                                   message="At least two files are required for",
                                   closure="\nSee dgeapy.py intersect --help for more information")
    if not args.names:
        raise MissingArgumentError("<NAME>", closure="\nSee dgeapy.py intersect --help for more information")
    elif len(args.names) != len(args.files):
        raise InvalidArgumentError(invalid_argument="-t <TABLE> -t <TABLE> [...] -n <NAME> -n <NAME> [...]",
                                   message="The number of names must match the number of files",
                                   closure="\nSee dgeapy.py intersect --help for more information")

    for f in args.files:
        if not os.path.isfile(f):
            raise FileNotFoundError(file_name=f)

    for frmt in args.formats:
        valid_formats = ["png", "svg", "pdf", "jpg", "jpeg"]
        if frmt not in valid_formats:
            raise InvalidArgumentError(invalid_argument=format,
                                       message=f"{TermMsg.ERROR}: Invalid format",
                                       closure=f"\nValid formats are: {Color.YELLOW}{', '.join(valid_formats)}{Color.RESET}")


def read_and_validate_tables(args: argparse.Namespace) -> list[pd.DataFrame]:
    """
    Try to read the input tables and perform the required validations:
        - Check if the table format is supported
        - Check if the index column contains NaN values
        - Check if the index column contains duplicated indexes

    Parameters:
        args (argparse.Namespace): The arguments.

    Returns:
        list[pd.DataFrame]: A list of pandas dataframes.

    Raises:
        UnsupportedFileFormatError: If the given table format is not supported
                                    or if the table could not be read.
        NanValuesInIndexColumnError: If the index column contains NaN values.
        DuplicatesInIndexColumnError: If the index column contains duplicated indexes.
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


    # Attempt to read each table
    tables = []
    for f in args.files:
        try:
            table = read_table(f, args.nan_values)
        except UnsupportedFileFormatError as e:
            raise e

        tables.append(table)


    # Look for NaN values in index columnn or duplicated indexes
    for idx, tbl in enumerate(tables):
        if not tbl[tbl[args.index_column].isna()].empty:
            raise NanValuesInIndexColumnError(file_index=idx)

        if not tbl[args.index_column].is_unique:
            duplicates = tbl[tbl[args.index_column].duplicated()][args.index_column].values.tolist()
            raise DuplicatesInIndexColumnError(file_index=idx,
                                               duplicated_indxs=duplicates)

        tbl.set_index(args.index_column, inplace=True)

    # Exclude indexes matching patterns
    if args.exclude:
        for pattern in args.exclude:
            for idx, tbl in enumerate(tables):
                has_pattern = tbl.index.str.contains(pattern)
                info_msg = (f"{TermMsg.INFO}: {Color.YELLOW}({has_pattern.sum()}{Color.RESET}) "
                            + f"indexes matching pattern '{Color.YELLOW}{pattern}{Color.RESET}' "
                            + "have been excluded")
                eprint(info_msg)

                # Actually excluding matching indexes
                tables[idx] = tbl[~has_pattern]


    return tables


def generate_output(args: argparse.Namespace, tables: list[pd.DataFrame]) -> None:

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

        Example usage:
            combinations = ["01", "10", "11"]
            dfs = [pd.DataFrame(index=[1, 2, 3]), pd.DataFrame(index=[2, 3, 4])]
            intersections = obtain_intersections(combinations, dfs)
            print(intersections)

            >> {'01': {4}, '10': {1}, '11': {2, 3}}
        """

        # Create a dictionary of indexes from each dataframe
        index_sets = [set(df.index) for df in dfs]

        intersections = {}

        for combination in combinations:
            included_sets = [index_sets[i] for i in range(len(combination)) if combination[i] == "1"]
            excluded_sets = [index_sets[i] for i in range(len(combination)) if combination[i] == "0"]

            if included_sets:
                intersection = set.intersection(*included_sets)
            else:
                intersection = set()

            for excluded_set in excluded_sets:
                intersection -= excluded_set

            if intersection:
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


    def generate_upset_plot(dfs_list_dict: dict[str, list[any]], output_directory: str, formats: list[str]) -> None:
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


    n_files = len(tables)

    # Will contain pd.DataFrames as values and names as keys
    # Will be used for generating intersections tables
    table_dict = {}
    # Will contain sets of indexes as values and names as keys
    # Will be used for generating venn diagrams
    venn_set_dict = {}
    # Will contain lists of indexes as values and names as keys
    # Will be used for generating upset plots
    upset_list_dict = {}
    for name, tbl in zip(args.names, tables):
        table_dict[name] = tbl
        venn_set_dict[name] = set(tbl.index.values.tolist())
        upset_list_dict[name] = tbl.index.values.tolist()

    intersections = obtain_intersctions(
                        combinations=obtain_all_possible_combinations(n_files),
                        dfs=tables)

    # Generating intersections dataframes
    intersections_dfs = generate_intersections_dfs(intersections, table_dict, args.index_column)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    match args.files[0].split(".")[-1]:
        case "csv":
            for k, v in intersections_dfs.items():
                v.to_csv(f"{args.output}/{k}.tsv", sep="\t")
        case "tsv":
            for k, v in intersections_dfs.items():
                v.to_csv(f"{args.output}/{k}.tsv", sep="\t")
        case "xlsx":
            for k, v in intersections_dfs.items():
                v.to_excel(f"{args.output}/{k}.xlsx")

    match n_files:
        case 2:
            generate_venn2_diagrams(venn_set_dict,
                                    args.output,
                                    args.formats)
            generate_upset_plot(upset_list_dict,
                                args.output,
                                args.formats)
        case 3:
            generate_venn3_diagrams(venn_set_dict,
                                    args.output,
                                    args.formats)
            generate_upset_plot(upset_list_dict,
                                args.output,
                                args.formats)
        case _:
            generate_upset_plot(upset_list_dict,
                                args.output,
                                args.formats)


def main():

    args = setup_parser().parse_args()

    # Validate arguments
    try:
        check_args(args)
    except (MissingArgumentError, InvalidArgumentError, FileNotFoundError) as e:
        eprint(f"{e}")
        sys.exit(1)

    try:
        tables = read_and_validate_tables(args)
    except (UnsupportedFileFormatError, NanValuesInIndexColumnError,
            DuplicatesInIndexColumnError) as e:
        eprint(f"{e}")
        sys.exit(1)

    generate_output(args, tables)

    eprint(f"{TermMsg.INFO}: Analysis done. See you again soon!")

    sys.exit(0)


if __name__ == "__main__":
    main()



