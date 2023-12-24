#!/usr/bin/env python3

"""
Script : dgeapy.py
Description : Differential Gene Expression Analyisis in Python at different levels.
Author : Joan Lluis Pons Ramon
Email : joanlluispons@gmail.com
"""


import os
import sys
import subprocess


__author__ = "Joan Lluis Pons Ramon"
__version__ = "0.0.0"


def main():

    description = """
Differential Gene Expression Analyisis in Python at different levels.

Usage: python dgeapy.py <COMMAND> [OPTIONS]

Commands:
    analyze           Perform differential gene expression analyisis
    intersect         Find intersections between indexes of multiple files

Options:
    -h, --help        Show this help message and exit
    -v, --version     Show version number and exit

Examples:
    python dgeapy.py analyze -h
    python dgeapy.py analyze --version

For more information, see <https://github.com/jllpons/dgeapy>.
    """

    arg_len = len(sys.argv)
    if arg_len == 1:
        print(description, file=sys.stderr)
        print("error: no command was given.", file=sys.stderr)
        sys.exit(1)

    elif arg_len > 1:
        cmd = sys.argv[1]
        dgeapy_path = os.path.dirname(os.path.realpath(__file__))

        if cmd == '-h' or cmd == '--help':
            print(description, file=sys.stderr)
            sys.exit(0)

        elif cmd == '-v' or cmd == '--version':
            print(f"dgeapy version {__version__}", file=sys.stderr)
            sys.exit(0)

        elif cmd == 'analyze':
            subcmd = ['python', f'{dgeapy_path}/dgeapy_analyze.py',] + sys.argv[2:]
            subprocess.run(subcmd)

        elif cmd == 'intersect':
            subcmd = ['python', f'{dgeapy_path}/dgeapy_intersect.py',] + sys.argv[2:]
            subprocess.run(subcmd)

        else:
            print(description, file=sys.stderr)
            print(f"error: command '{cmd}' not found.", file=sys.stderr)
            sys.exit(1)


if __name__ == '__main__':
    main()
