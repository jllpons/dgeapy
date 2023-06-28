# dgeapy: Differential Gene Expression Analysis in Python

dgeapy.py is a set of tools that try to analyze RNAseq data at different levels.

Available scripts are:

- [dgeapy dgea](###-dgeapy-dgea) 

## Dependencies

Python (>= 3.9)

```shell
conda create --name dgeapy python=3.9
```

- Data analysis:
    - [pandas](<https://pypi.org/project/pandas/>): dataframe analysis and
      manipulation
    - [numpy](<https://pypi.org/project/numpy/>): computing
    - [openpyxl](<https://pypi.org/project/openpyxl/>): engine for reading/writting xlsx files
- Data visualization:
    - [matplotlib](<https://pypi.org/project/matplotlib/>): low level manipulations.
    - [seaborn](<https://pypi.org/project/seaborn/>): high level manipulations.

```shell
pip install pandas numpy openpyxl matplotlib seaborn
```
## Usage

```
./dgeapy.py -h

dgeapy: Differential Gene Expression Analyisis in Python at different levels.

usage: dgeapy.py <command> [options]

    dgea    differential gene expression analyisis

    options:
        -h, --help
```

### dgeapy.py dgea

```
./dgeapy.py dgea -h
usage: dgeapy.py dgea TABLE [options]

Perform Differential Gene Expression Analysis (DGEA) by determining the differentially expressed genes from a dataframe. It takes as input a table in CSV, TSV, or XLSX
format containing gene expression data. The script applies thresholds for adjusted p-values and fold changes to identify significant gene expression changes. Generates bar
plots and volcano plots to visualize the results. The output includes the modified dataframe with added columns for fold change and gene regulation, as well as the generated
plots saved in the specified output directory.

positional arguments:
  TABLE                 path to the CSV, TSV or XLSX file

optional arguments:
  -h, --help            show this help message and exit
  -o STR, -output-directory STR
                        output directory [Default: $CWD/dgeapy_output]
  --padj FLOAT          adjusted p-value threshold. LESS THAN THRESHOLD WILL BE APPLIED [Default: 0.05]
  --fc FLOAT            fold change threshold. ABSOLUTE VALUE EQUAL OR MORE THAN THRESHOLD WILL BE APPLIED [Default: 1.5]
  --formats [STR]       plot formats [Default: ['png']]
  --exclude [STR]       string pattern to match against indexes. Matched indexes are excluded
  --nan-values [STR]    strings to recognize as NaN values. Transcripts with NaN padj or NaN fold change will be excluded [Default: ["", "--"]]
  --keep-duplicated     if passed, keep duplicate index values [Default: False]
  --index-column STR    name of the column that will be used as index [Default: index]
  --log2fc-column STR   name of the column containing the log2 Fold Change values [Default: log2_fold_change]
  --padj-column STR     name of the column containing the p-adjusted values [Default: padj]
```

