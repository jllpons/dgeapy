#!/usr/bin/env python3

"""
Given a list of dataframes, compute all possible intersections of the indexes
of the dataframes. Save one dataframe for each intersection and generate venn
diagrams and upset plots for data visualization.
"""


import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_unweighted, venn3, venn3_unweighted
from upsetplot import from_contents, UpSet
import pandas as pd


# Defalut colum names
defalut_index_column_name = "index"

# Defalut plot formats
defalut_plot_formats = ["png", "pdf"]


def obtain_all_possible_combinations(n):
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


def obtain_intersctions(combinations, dfs):
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


def generate_intersections_dfs(intersections, dfs_dict, index_column):
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


def generate_venn2_diagrams(dfs_dict, output_directory, formats):
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


def generate_venn3_diagrams(dfs_set_dict, output_directory, formats):
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


def generate_upset_plot(dfs_list_dict, output_directory, formats):
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

    description = """
    Given a list of data files, compute all the possible intersections between them.
    Then, generate two tables (TSV and XLSX) for each intersection. If the number of
    data files is less than or equal to 3, generate Venn diagrams and an UpSet plot.
    If the number of data files is greater than 3, generate only an UpSet plot.
    """

    parser = argparse.ArgumentParser(
            description=description,
            usage="dgeapy.py intersections <file1> <file2> <file3> ... name_1 name_2 name_3 ... [options]",
            )

    parser.add_argument(
            "-f", "--files",
            metavar="<file_1> <file_2> ...",
            nargs="+",
            type=str,
            help="data files to be processed",
            )
    parser.add_argument(
            "-n", "--names",
            metavar="<name_1> <name_2> ...",
            nargs="+",
            type=str,
            help="names of the data files that will be used in the plots and tables",
            )
    parser.add_argument(
            "-o", "--output_directory",
            metavar="PATH",
            default=f"{os.getcwd()}/dgeapy_intersections_output",
            type=str,
            help="output directory [Default: $CWD/dgeapy_intersections_output]",
            )
    parser.add_argument(
            "-i", "--index_column",
            metavar="STR",
            default=defalut_index_column_name,
            type=str,
            help="name of the index column",
            )
    parser.add_argument(
            "--formats",
            metavar="STR",
            nargs="+",
            default=defalut_plot_formats,
            type=str,
            action="extend",
            help=f"output formats for the plots [Default: {defalut_plot_formats}]",
            )
    parser.add_argument(
            "--nan-values",
            metavar="STR",
            nargs="+",
            default=["", "--", "NA"],
            action="extend",
            type=str,
            help="strings to recognize as NaN values in index column Default: ['', '--', 'NA']",
            )
    parser.add_argument(
            "--exclude",
            metavar="STR",
            nargs="+",
            default=[],
            action="extend",
            type=str,
            help="string patterns to exclude from the index column",
            )

    args = parser.parse_args()

    if not args.files:
        parser.print_help()
        sys.exit("\n** ERROR: No files were provided **")

    if not args.names:
        parser.print_help()
        sys.exit("\n** ERROR: No names were provided **")

    files = args.files
    names = args.names
    output_directory = args.output_directory

    index_column = args.index_column
    nan_values = args.nan_values
    pattern_exclude = args.exclude

    formats = args.formats

    dfs = []
    for f in files:
        if not os.path.isfile(f):
            sys.exit(f"** ERROR: could not find {f} **")

        try:

            if f.endswith(".csv"):
                df = pd.read_csv(f, na_values=nan_values)
                dfs.append(df)
            elif f.endswith(".tsv"):
                df = pd.read_csv(f, sep="\t", na_values=nan_values)
                dfs.append(df)
            elif f.endswith(".xlsx"):
                df = pd.read_excel(f,  na_values=nan_values)
                dfs.append(df)
            else:
                sys.exit(f"** ERROR: unsupported format for {f}. Supported formats "
                         + "are CSV, TSV and XLSX **")
        except:
            sys.exit(f"** ERROR: could not read {f} **")

    for i, df in enumerate(dfs):
        if index_column not in df.columns.values.tolist():
            print(df.columns.values.tolist())
            sys.exit(f"** ERROR: could not find {index_column} in file_{i} **")
        if not df[index_column].is_unique:
            sys.exit(f"** ERROR: {index_column} contains duplicated values in file_{i} **")
        if not df[df[index_column].isna()].empty:
            sys.exit(f"** ERROR: {index_column} contains NaN values in file_{i} **")

        df.set_index(index_column, inplace=True)

    # Exclude patterns from index
    if pattern_exclude:
        new_dfs = []
        for i, df in enumerate(dfs):
            for pattern in pattern_exclude:
                df = df[~df.index.str.contains(pattern)]
            new_dfs.append(df)
        dfs = new_dfs


    n_files = len(dfs)

    intersections = obtain_intersctions(
                        combinations=obtain_all_possible_combinations(n_files),
                        dfs=dfs,
                        )

    # Dictionary of dataframes (for tables)
    dfs_dict = {}
    # Dictionary of sets of indexes (for venn diagram functions)
    dfs_set_dict = {}
    # Dictionary of lists of indexes (for upset plot function)
    dfs_list_dict = {}
    for name, df in zip(names, dfs):
        dfs_dict[name] = df
        dfs_set_dict[name] = set(df.index.values.tolist())
        dfs_list_dict[name] = df.index.values.tolist()

    # Generating intersections dataframes
    intersections_dfs = generate_intersections_dfs(intersections, dfs_dict, index_column)
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)
    for key, value in intersections_dfs.items():
        value.to_csv(f"{output_directory}/{key}.tsv", sep="\t")
        value.to_excel(f"{output_directory}/{key}.xlsx")

    fig_output_directory = f"{output_directory}/fig"
    if not os.path.isdir(fig_output_directory):
        os.mkdir(fig_output_directory)

    if n_files == 2:
        generate_venn2_diagrams(
                dfs_set_dict,
                fig_output_directory,
                formats,
                )
        generate_upset_plot(
                dfs_list_dict,
                fig_output_directory,
                formats,
                )

    elif n_files == 3:
        generate_venn3_diagrams(
                dfs_set_dict,
                fig_output_directory,
                formats,
                )
        generate_upset_plot(
                dfs_list_dict,
                fig_output_directory,
                formats,
                )

    else:
        generate_upset_plot(
                dfs_list_dict,
                fig_output_directory,
                formats,
                )


if __name__ == "__main__":
    main()



