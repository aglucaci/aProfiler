#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# Imports
# =============================================================================
#from AlignmentProfiler import *
#import AlignmentProfiler as AP
from AlignmentProfiler import MSA_Process
import argparse
import os
import sys
import pathlib 
from tabulate import tabulate

# =============================================================================
# Declares
# =============================================================================
parser = argparse.ArgumentParser(description='Process command line arguments')

parser.add_argument('--alignment', metavar='A', type=str, required=True,
                    help='Input multiple sequence alignment')
                    
parser.add_argument('--output', metavar='O', type=str, required=False,
                    help='Specify a location for the output CSV')

args = parser.parse_args()

#OutputDirectory = os.path.join("..", "..", "tables")
OutputDirectory = os.path.join("tables")

if args.output == None:
    OutputCSV = os.path.join(OutputDirectory,
                             os.path.basename(args.alignment) + "_output.csv")
else:
    OutputCSV = args.output
#end if

print("# Alignment file is:", args.alignment)
print("# Output will be saved to:", OutputCSV)

# =============================================================================
# Main
# =============================================================================
print()
print("#", "".join(["="]*78))
print("# Starting to profile multiple sequence alignment:", args.alignment)
print("#", "".join(["="]*78))
print()

if os.path.exists(args.alignment) and os.path.getsize(args.alignment) > 0:
    _MSA = MSA_Process(args.alignment)
    print()
    print("# Loading complete on an alignment with", _MSA.NumSequences, "sequences, and", _MSA.NumCodonSites, "codon sites")
    print()
else:
    print("# ERROR: the file is either empty or does not exist")
#end if


# =============================================================================
# Report to user
# =============================================================================

print("#", "".join(["="]*78))
print("# Reporting alignment statistics")
print("#", "".join(["="]*78))
print()

# Invariant Sites
print("# The alignment has", _MSA.NumInvariantCodonSites, "invariant sites", "(", (_MSA.NumInvariantCodonSites / _MSA.NumCodonSites) * 100, "% ) ")

# Gaps
#print("# The alignment has an average gappiness of", 0, "with a range of", 0, "to", 0)

# Nucleotide Frequencies
print("# The alignment the following nucleotide frequencies...")
print()
print(tabulate([
    ['Adenine (A)',  _MSA.NTFrequencies["NucleotideFrequencies"]["A"]],
    ['Thymine (T)',  _MSA.NTFrequencies["NucleotideFrequencies"]["T"]],
    ['Guanine (G)',  _MSA.NTFrequencies["NucleotideFrequencies"]["G"]],
    ['Cytosine (C)', _MSA.NTFrequencies["NucleotideFrequencies"]["C"]],
    ['Any Nucleotide (N)', _MSA.NTFrequencies["NucleotideFrequencies"]["N"]]],
    headers=['Nucleotide', 'Frequency']))


# =============================================================================
# End of file
# =============================================================================
