#!/usr/bin/env python3

"""
dgeapy: a script that tries to analyse Differential Gene Expression (DGE)
data at a different levels.
"""

import os
import sys
import subprocess


def main():

    description = """
dgeapy: Differential Gene Expression Analyisis at different levels.

usage: dgeapy.py <command> [options]

    dgea    differential gene expression analyisis

    options:
        -h, --help
    """

    arg_len = len(sys.argv)
    if arg_len == 1:
        print(description)

    elif arg_len > 1:
        cmd = sys.argv[1]
        dgeapy_path = os.path.dirname(os.path.realpath(__file__))

        if cmd == '-h' or cmd == '--help':
            print(description)

        elif cmd == 'dgea':
            subcmd = ['python', f'{dgeapy_path}/dgeapy_dgea.py',] + sys.argv[2:]
            subprocess.run(subcmd)

        else:
            print(description)


if __name__ == '__main__':
    main()
