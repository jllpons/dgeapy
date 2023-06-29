# dgeapy: Differential Gene Expression Analysis in Python

dgeapy.py is a set of tools that try to analyze RNAseq data at different levels.

Available scripts are:

- [dgeapy dgea](#dgeapy-dgea) 

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

Determine the differentially expressed genes from a dataframe.

1. Takes a table in CSV, TSV or XLSX format as input.
2. Checks for any NaN values in index, p-adjusted and log2 fold change columns.
3. Indexes each row of the table according to `--index-colum`. 
4. Checks for duplicated values in index.
5. Excludes indexes that contain specific patterns provided with the `--exclude` argv.
6. Adds fold change and regulation columns.
7. **Determines the differentially expressed genes** (DEG). A gene is considered as DEG when it's p-adjusted value is less than the provided threshold in `--padj` and it's fold change absoulte value is equal or greater than the provided threshold in `--fc`.
8. Creates an output directory `-o` where it stores a table containing the DEG `deg.table`, a table containing only the upregulated genes `upregulated.table` and another table containing the downregulated genes `downregulated.table` in both TSV and XLSX formats.
9. Creates a `fig` directory inside the output directory containing a generated bar plot (Fig. 1)
   and volcano plot (Fig. 2).

```
./dgeapy.py dgea -h
usage: dgeapy.py dgea TABLE [options]

Perform Differential Gene Expression Analysis (DGEA) by determining the differentially expressed
genes from a dataframe. It takes as input a table in CSV, TSV, or XLSX format containing gene
expression data. The script applies thresholds for adjusted p-values and fold changes to
identify significant gene expression changes. Generates bar plots and volcano plots to visualize
the results. The output includes the modified dataframe with added columns for fold change and
gene regulation, as well as the generated plots saved in the specified output directory.

positional arguments:
  TABLE                 path to the CSV, TSV or XLSX file

optional arguments:
  -h, --help            show this help message and exit
  -o STR, -output-directory STR
                        output directory [Default: $CWD/dgeapy_output]
  --padj FLOAT          adjusted p-value threshold. LESS THAN THRESHOLD WILL BE APPLIED [Default: 0.05]
  --fc FLOAT            fold change threshold. ABSOLUTE VALUE EQUAL OR MORE THAN THRESHOLD WILL BE APPLIED [Default: 1.5]
  --formats STR       plot formats [Default: ['png']]
  --exclude STR       string pattern to match against indexes. Matched indexes are excluded
  --nan-values STR    strings to recognize as NaN values. Transcripts with NaN padj or NaN fold change will be excluded [Default: ["", "--"]]
  --keep-duplicated     if passed, keep duplicate index values [Default: False]
  --index-column STR    name of the column that will be used as index [Default: index]
  --log2fc-column STR   name of the column containing the log2 Fold Change values [Default: log2_fold_change]
  --padj-column STR     name of the column containing the p-adjusted values [Default: padj]
```

#### Example usage:

```shell
./dgeapy/dgeapy.py dgea example/data/dgeapy_dgea_example.xlsx -o example/dgeapy_output --log2fc-column log2FoldChange --nan-values NA
```

Tables and figures in `example/dgeapy_output`.

![**Figure 1**: Bar plot generated using the exmple data.](data/dgeapy_output/fig/barplot.png)

![**Figure 2**: Volcano plot generated using the exmple data.](data/dgeapy_output/fig/volcano.png)

Example data can be downloaded at [GSE206442](<https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE206442&format=file&file=GSE206442%5FGIBERT%5F01%5Fnew%5Fannot%5Fwo%5Foutlier%5FSTAT%5Fvs%5FLOG%5Fresults%2Exlsx>).

