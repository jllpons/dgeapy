# dgeapy: Differential Gene Expression Analysis in Python

DGEAPY is a Python toolkit for analyzing RNAseq data, focusing on differential gene expression and intersections among datasets.

## Installation

### Setting up the Environment

Create a dedicated environment for DGEAPY using Conda:

```shell
conda create --name dgeapy python=3.10
conda activate dgeapy
```

### Installing Dependencies

Install all necessary libraries with Conda:

```shell
conda install pandas numpy openpyxl matplotlib seaborn matplotlib-venn UpSetPlot
```

## Usage

Available scripts are:

- [dgeapy.py analyze](#dgeapy-analyze) 
- [dgeapy.py intersect](#dgeapy-intersect)

```
python dgeapy/dgeapy.py -h

Differential Gene Expression Analyisis in Python at different levels.

Usage: python dgeapy.py <COMMAND> [OPTIONS]

Commands:
    analyze           Perform differential gene expression analyisis
    intersect         Find intersections between indexes of multiple files

Options:
    -h, --help        Show this help message and exit
    -v, --version     Show version number and exit
```

### dgeapy analyze

Determine the differentially expressed genes from a dataframe.

1. Input a table in CSV, TSV, or XLSX format.
2. Verify and clean the data by checking for NaN values, duplicated values in the index, and excluding indexes with specific patterns using `--exclude`.
3. Utilize `--index-column` to index each row and add fold change and regulation columns.
4. **Identify differentially expressed genes (DEG)** by applying thresholds for p-adjusted value (`--padj`) and fold change absolute value (`--fc`).
5. Output three table (DEGs, upregulated and downregulated) and two figures (a bar plot and a volcano plot).

```
python dgeapy/dgeapy.py analyze -h
usage: dgeapy.py analyze <TABLE> [OPTIONS]

Differential Gene Expression Analysis.
Generates tables with differentially expressed plots to visualize the results.


positional arguments:
  <TABLE>                  Path to the gene expression data file (CSV, TSV, XLSX).

options:
  -h, --help               show this help message and exit
  -o, --output DIR         Specify the output directory (default: cwd).
  -p, --padj FLOAT         Adjusted p-value threshold for significance (default: 0.05).
  -f, --fold-change FLOAT  Fold change threshold for significance (default: 1.5).
  -F, --formats [STR]      Output formats for plots (e.g. svg, pdf) (default: ['png']).
  -e, --exclude [STR]      Exclude indexes matching specified patterns.
  -N, --nan-values [STR]   Strings to recognize as NaN (default: ['', '--']).
  -k, --keep-duplicated    Keep duplicated indexes (default: False).
  -I, --index-column STR   Column name for index (default: index).
  -L, --log2fc-column STR  Column name for log2 Fold Change (default: log2_fold_change).
  -P, --p-column STR       Column name for adjusted p-values (default: padj).
```

#### Usage example:

```shell
python dgeapy/dgeapy.py analyze example/data/dgeapy_dgea_example.xlsx -o example/analyze_output -L log2FoldChange -N NA
```

Output tables and figures can be found in `example/analyze_output`.

![**Figure 1**: Bar plot generated using the exmple data.](example/analyze_output/barplot.png)

![**Figure 2**: Volcano plot generated using the exmple data.](example/analyze_output/volcano.png)

Example data can be downloaded at [GSE206442](<https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE206442&format=file&file=GSE206442%5FGIBERT%5F01%5Fnew%5Fannot%5Fwo%5Foutlier%5FSTAT%5Fvs%5FLOG%5Fresults%2Exlsx>).

### dgeapy intersect

Compute intersections of indexes among a list of dataframes.

1. Takes multiple dataframes and their assigned name.
2. Checks for any NaN values, null values or duplicaded values in the indexes.
3. Excludes indexes that contain specific patterns provided with the `--exclude` argv.
4. **Computes all of the possible intersections between the indexes of the provided dataframes.**
5. Generates a TSV and a XLSX file for each non-empty intersection.
6. Generates a weighted and an unweighted Venn Diagram if the n of provided
   tables is 3 or less.
7. Generates an UpSet Plot for a better visualitzation of the present and
   missing intersections.

```
./dgeapy/dgeapy.py intersections -h
usage: dgeapy.py intersections <file1> <file2> <file3> ... name_1 name_2 name_3 ... [options]

Given a list of data files, compute all the possible intersections between them.
Then, generate two tables (TSV and XLSX) for each intersection. If the number of 
data files is less than or equal to 3, generate Venn diagrams and an UpSet plot.
If the number of data files is greater than 3, generate only an UpSet plot.

optional arguments:
  -h, --help            show this help message and exit
  -f <file_1> <file_2> ... [<file_1> <file_2> ... ...], --files <file_1> <file_2> ... [<file_1> <file_2> ... ...]
                        data files to be processed
  -n <name_1> <name_2> ... [<name_1> <name_2> ... ...], --names <name_1> <name_2> ... [<name_1> <name_2> ... ...]
                        names of the data files that will be used in the plots and tables
  -o PATH, --output_directory PATH
                        output directory [Default: $CWD/dgeapy_intersections_output]
  -i STR, --index_column STR
                        name of the index column
  --formats STR [STR ...]
                        output formats for the plots [Default: ['png', 'pdf']]
  --nan-values STR [STR ...]
                        strings to recognize as NaN values in index column Default: ['', '--', 'NA']
  --exclude STR [STR ...]
                        string patterns to exclude from the index column
```

#### Usage example:

Intersection analysis between 3 files:

```shell
./dgeapy/dgeapy.py intersections -f test/dgeapy_intersections_test/data/set1-20.xlsx test/dgeapy_intersections_test/data/set4-23.xlsx test/dgeapy_intersections_test/data/set7-26.xlsx -n set-1-20 set-4-23 set-7-26 -o test/dgeapy_intersections_test/3_sets_test
```
Results can be found in `test/dgeapy_intersections_test/3_sets_test`

Example of a generated Venn Diagram: 

![**Figure 3**: Venn diagram using example
data](test/dgeapy_intersections_test/3_sets_test/fig/venn3_unweighted.png)

Intersection analysis between 4 files:

Results can be found in `test/dgeapy_intersections_test/4_sets_test`

Example of a generated UpSet Plot:

![**Figure 4**: Venn diagram using example
data](test/dgeapy_intersections_test/4_sets_test/fig/upset.png)

```shell
./dgeapy/dgeapy.py intersections -f test/dgeapy_intersections_test/data/set1-20.xlsx test/dgeapy_intersections_test/data/set4-23.xlsx test/dgeapy_intersections_test/data/set7-26.xlsx test/dgeapy_intersections_test/data/set16-35.xlsx -n set-1-20 set-4-23 set-7-26 set-16-35 -o test/dgeapy_intersections_test/4_sets_test
```
