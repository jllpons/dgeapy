#!/usr/bin/env python3

"""
Different utilities for dgeapy.
"""

import argparse

__author__ = "Joan Lluis Pons Ramon"
__email__ = "joanlluis@gmail.com"
__version__ = "0.0.0"

class Color:
    """
    Color class for terminal output.
    """

    GREEN = "\033[32m"
    RED = "\033[31m"
    YELLOW = "\033[33m"

    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"

    RESET = "\033[0m"

class TermMsg:
    """
    Templates for terminal messages.
    """

    ERROR = f"[dgeapy] {Color.RED}Error{Color.RESET}"
    WARNING = f"[dgeapy] {Color.RED}Warning{Color.RESET}"
    INFO = f"[dgeapy] {Color.GREEN}Info{Color.RESET}"


# <https://stackoverflow.com/questions/18275023/dont-show-long-options-twice-in-print-help-from-argparse>
# I have also added `argparse.RawTextHelpFormatter` to the inheritance
# because I want to be able to use "\n" in the description and epilog.
class CustomHelpFormatter(argparse.RawTextHelpFormatter, argparse.HelpFormatter):
    def __init__(self, prog):
        # Initialize with super from RawTextHelpFormatter and also set our custom widths
        super().__init__(prog, max_help_position=40)

    def _format_action_invocation(self, action):
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string



