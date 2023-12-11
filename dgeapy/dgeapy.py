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
Differential Gene Expression Analyisis in Python at different levels.

Usage: python dgeapy.py <COMMAND> [OPTIONS]

Commands:
    dgea            differential gene expression analyisis
    intersections   find intersections between the indexes of n files

Options:
    -h, --help      show this help message and exit
    """

    arg_len = len(sys.argv)
    if arg_len == 1:
        print(description, file=sys.stderr)
        print('error: missing subcommand', file=sys.stderr)
        sys.exit(1)

    elif arg_len > 1:
        cmd = sys.argv[1]
        dgeapy_path = os.path.dirname(os.path.realpath(__file__))

        if cmd == '-h' or cmd == '--help':
            print(description, file=sys.stderr)

        elif cmd == 'dgea':
            subcmd = ['python', f'{dgeapy_path}/dgeapy_dgea.py',] + sys.argv[2:]
            subprocess.run(subcmd)

        elif cmd == 'intersections':
            subcmd = ['python', f'{dgeapy_path}/dgeapy_intersections.py',] + sys.argv[2:]
            subprocess.run(subcmd)

        else:
            print(description, file=sys.stderr)
            print(f'error: unknown subcommand "{cmd}"', file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
