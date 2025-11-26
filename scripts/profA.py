#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""


"""

# =============================================================================
# Imports
# =============================================================================

import argparse
import os
import sys
import pathlib
#from tabulate import tabulate

# =============================================================================
# Parse command line arguments
# =============================================================================
parser = argparse.ArgumentParser(description='')

parser.add_argument('--input', metavar='I', type=str, required=True,
                    help='Input multiple sequence alignment')
parser.add_argument('--output', metavar='O', type=str, required=False,
                    help='Specify a directory for outputs')
parser.add_argument('--file-type', metavar='T', type=str, required=False,
                    help='Specify a directory for outputs')
args = parser.parse_args()

# =============================================================================
# Declares
# =============================================================================

_outputDirectory = os.path.join("tables")

if args.output == None:
    outputCSV = os.path.join(_outputDirectory,
                             os.path.basename(args.input) + "_output.csv")
else:
    outputCSV = args.output
#end if

print("# Alignment file is:", args.input)
print("# Output will be saved to:", outputCSV)

# =============================================================================
# Main
# =============================================================================
#print()
#print("#", "".join(["="]*78))
print("# [INFO] Starting to profile multiple sequence alignment:", args.input)
print("# [INFO]", args.input)
#print("#", "".join(["="]*78))
#print()



# =============================================================================
# End of file
# =============================================================================
