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

from dgeapy_utils import ( __version__, Color, TermMsg, eprint)


def main():

    description = """
Differential Gene Expression Analyisis in Python at different levels.

Usage: python dgeapy.py <COMMAND> [OPTIONS]

Commands:
    analyze           Perform differential gene expression analyisis
    intersect         Find intersections between indexes of multiple files

Options:
    -h, --help        Show this help message and exit
    -V, --version     Show version number and exit

Examples:
    python dgeapy.py analyze -h
    python dgeapy.py analyze <mygenes.csv> -f 2.0 -p 0.01
    python dgeapy.py intersect -f <mutA.csv> -f <mutB.csv> -n "Mutant A" -n "Mutant B" -i gene_id

For more information, see <https://github.com/jllpons/dgeapy>.
    """

    arg_len = len(sys.argv)
    if arg_len == 1:
        eprint(description)
        eprint(f"{TermMsg.ERROR}: no command was given.")
        sys.exit(1)

    elif arg_len > 1:
        cmd = sys.argv[1]
        dgeapy_path = os.path.dirname(os.path.realpath(__file__))

        if cmd == '-h' or cmd == '--help':
            print(description)
            sys.exit(0)

        elif cmd == '-V' or cmd == '--version':
            print(f"dgeapy version {Color.YELLOW}{__version__}{Color.RESET}")
            sys.exit(0)

        elif cmd == 'analyze':
            subcmd = ['python', f'{dgeapy_path}/dgeapy_analyze.py',] + sys.argv[2:]
            subprocess.run(subcmd)

        elif cmd == 'intersect':
            subcmd = ['python', f'{dgeapy_path}/dgeapy_intersect.py',] + sys.argv[2:]
            subprocess.run(subcmd)

        else:
            eprint(description)
            eprint(f"{TermMsg.ERROR}: unrecognized command '{Color.YELLOW}{cmd}{Color.RESET}'.")
            sys.exit(1)


if __name__ == '__main__':
    main()
