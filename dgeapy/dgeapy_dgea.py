#!/usr/bin/env python3

"""
Perform Differential Gene Expression Analysis (DGEA) by determining the
differentially expressed genes from a dataframe. It takes as input a table
in CSV, TSV, or XLSX format containing gene expression data. The script
applies thresholds for adjusted p-values and fold changes to identify
significant gene expression changes.

Generates bar plots and volcano plots to visualize the results.

The output includes the modified dataframe with added columns for fold
change and gene regulation, as well as the generated plots saved in the
specified output directory.
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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

    column_names.insert(fold_change_column_position, fold_change_column_name)

    return df.reindex(columns=column_names)


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

    df.loc[df[log2_fold_change_column_name] > 0, regulation_column_name] = "Upregulated"
    df.loc[df[log2_fold_change_column_name] < 0, regulation_column_name] = "Downregulated"


    column_names.insert(regulation_column_position, regulation_column_name)

    return df.reindex(columns=column_names)


def mk_bar_plot(data, output_directory, plot_formats):
    """
    Generate a bar plot with count and percentage annotations.

    Parameters:
        data (pandas.DataFrame): The data for creating the bar plot.
        output_directory (str): The directory to save the generated plot.
        plot_formats (list): A list of plot formats to save the plot in.

    Returns:
        None
    """

    plt.figure(figsize=(11,6))

    count_plot = sns.countplot(
                        data=data,
                        y="significance",
                        palette=("silver", "cornflowerblue", "indianred"),
                        order=["No significant", "Downregulated", "Upregulated"],
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

    fig = count_plot.get_figure()

    for format in plot_formats:
        plot_name = f"{output_directory}/barplot.{format}"
        fig.savefig(plot_name, format=format, dpi=300)

        if format == "png":
            plot_name = f"{output_directory}/barplot_transparent-bg.{format}"
            fig.savefig(plot_name, format=format, dpi=300, transparent=True)

    plt.close()


def mk_volcano_plot(
        data,
        log2_fold_change_column_name,
        padj_column_name,
        foldchange_threshold,
        padj_threshold,
        output_directory,
        plot_formats
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

    description = """
    Perform Differential Gene Expression Analysis (DGEA) by determining the
    differentially expressed genes from a dataframe. It takes as input a table
    in CSV, TSV, or XLSX format containing gene expression data. The script
    applies thresholds for adjusted p-values and fold changes to identify
    significant gene expression changes.

    Generates bar plots and volcano plots to visualize the results.\n

    The output includes the modified dataframe with added columns for fold
    change and gene regulation, as well as the generated plots saved in the
    specified output directory.
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

    upregulated_genes_df = diff_expressed_genes_df[diff_expressed_genes_df[regulation_column_name] == "Upregulated"]
    downregulated_genes_df = diff_expressed_genes_df[diff_expressed_genes_df[regulation_column_name] == "Downregulated"]

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    diff_expressed_genes_df.to_csv(f"{output_directory}/deg.tsv", sep="\t")
    upregulated_genes_df.to_csv(f"{output_directory}/upregulated.tsv", sep="\t")
    downregulated_genes_df.to_csv(f"{output_directory}/downregulated.tsv", sep="\t")

    diff_expressed_genes_df.to_excel(f"{output_directory}/deg.xlsx")
    upregulated_genes_df.to_excel(f"{output_directory}/upregulated.xlsx")
    downregulated_genes_df.to_excel(f"{output_directory}/downregulated.xlsx")

    # Preparing values for data representation
    df["significance"] = np.where(df.index.isin(diff_expressed_genes_df.index), "Significant", "No significant")
    df.loc[df["significance"] == "Significant", "significance"] = df[regulation_column_name]

    fig_directory = output_directory + "/fig"
    if not os.path.isdir(fig_directory):
        os.mkdir(fig_directory)

    mk_bar_plot(df, fig_directory, plot_formats)
    mk_volcano_plot(
            data=df,
            log2_fold_change_column_name=log2_fold_change_column_name,
            padj_column_name=padj_column_name,
            foldchange_threshold=fold_change_threshold,
            padj_threshold=padj_threshold,
            output_directory=fig_directory,
            plot_formats=plot_formats,
            )

if __name__ == "__main__":
    main()
