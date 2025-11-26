#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Author: Alexander G. Lucaci

@Description: AlignmentProfiler provides quick summary statistics
                  for a multiple sequence alignment.

"""

# =============================================================================
# Imports
# =============================================================================

import pandas as pd
#from pandas.plotting import scatter_matrix
#import random
#from matplotlib import pyplot as plt
import itertools
import numpy as np
#from matplotlib.pyplot import figure
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#import numpy as np
#from pandas import DataFrame
#import seaborn as sns
from tqdm import tqdm
import argparse
import os
import sys
import pathlib 
#from tabulate import tabulate

# =============================================================================
# Declares
# =============================================================================
"""
parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--input', metavar='I', type=str, required=True,
                    help='Input multiple sequence alignment')
parser.add_argument('--output', metavar='O', type=str, required=False,
                    help='Specify a location for the output csv')


args = parser.parse_args()

OutputDirectory = os.path.join("..", "..", "tables")

if args.output == None:
    #print([args.output])
    OutputCSV = os.path.join(OutputDirectory, 
                             os.path.basename(args.input) + "_output.csv")
else:
    OutputCSV = args.output
#end if
"""

#print(args.input)
#print(OutputCSV)

Universal_codon_table = {"Glycine": ["GGT", "GGC", "GGA", "GGG"], 
                         "Alanine": ["GCT", "GCC", "GCA", "GCG"],
                         "Valine": ["GTT", "GTC", "GTA", "GTG"],
                         "Leucine": ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
                         "Methionine": ["ATG"],
                         "Isoleucine":  ["ATT", "ATC", "ATA"],
                         "Phenylalanine": ["TTT", "TTC"] ,
                         "Tyrosine": ["TAT", "TAC"],
                         "Tryptophan": ["TGG"],
                         "Serine": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                         "Threonine": ["ACT", "ACC", "ACA", "ACG"],
                         "Cysteine": ["TGT", "TGC"],
                         "Proline": ["CCT", "CCC", "CCA", "CCG"],
                         "Asparagine": ["AAT", "AAC"],
                         "Glutamine":  ["CAA", "CAG"], 
                         "Lysine": ["AAA", "AAG"],
                         "Arginine": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], 
                         "Histidine": ["CAT", "CAC"],
                         "Aspartate": ["GAT", "GAC"], 
                         "Glutamate": ["GAA", "GAG"],
                         "Stop":  ["TAA", "TAG", "TGA"]}


ResultsDict = {}

# =============================================================================
# Classes
# =============================================================================

class Alignment():
  def __init__(self, _Alignment, MoleculeType="DNA"):
    File_Extension = pathlib.Path(_Alignment).suffix
    if File_Extension == ".nex":
        File_Type = "nexus"
    else:
        File_Type = "fasta"
    #end if
    
    #self.records = self.process_sites(FASTA_FILE)
    #self.CodonTable = CodonTable
    
    self.NumCodonSites = self.Get_NumCodonSites(_Alignment, File_Type)
    self.NumSequences = self.Get_NumSequences(_Alignment, File_Type)
    
    #self.NumSequences = self.Get_NumSequences()
    #self.alignment_type = "DNA" # Default
    
    self.MoleculeType = MoleculeType
    self.GapCharacter = "-"                    # Default
    self.DNA_Characters = ["T", "C", "G", "A"] # Default
    
    # Statistics
    self.NumInvariantCodonSites = self.Get_NumInvariantCodonSites(_Alignment, File_Type)
    
    self.Gappiness = {}
    
    self.NTFrequencies = self.Get_NTFrequencies(_Alignment, File_Type)
  #end method

  def Get_NumCodonSites(self, _Alignment, File_Type):
      sites = 0
      with open(_Alignment, "r") as handle:
          for record in SeqIO.parse(handle, File_Type):
              if len(str(record.seq)) % 3 == 0:
                  sites = int(len(str(record.seq)) / 3)
                  return sites
              else:
                  print("# ERROR: Number of sites is not a multiple of 3")
            #end if
          #end for
          
      #end with
      #return sites, sequences  
  #end method

  def Get_NumSequences(self, _Alignment, File_Type):
      sequences = 0
      with open(_Alignment, "r") as handle:
          for record in SeqIO.parse(handle, File_Type):
              sequences += 1
          #end for
      #end with
      return sequences  
  #end method
  
  def Get_NumInvariantCodonSites(self, _Alignment, File_Type):
      InvariantSites = 0
      #VariableSites = 0
      #loop over every site
      #for i in range(0, self.NumSites):
      print("# Checking alignment for invariant sites")
      for i in tqdm(range(self.NumCodonSites)):
            #print("# Checking site:", i + 1)
            column_data = []
            Variable = False
            with open(_Alignment, "r") as handle:
                for record in SeqIO.parse(handle, File_Type):
                    codon = record.seq[i * 3 : (i * 3) + 3]
                    #print(i + 1, codon)
                    if column_data == []:
                        column_data.append(codon)
                    else:
                        # column data is not empty
                        # check it
                        if column_data[0] != codon:
                            # Variable site
                            # Go to next site
                            #print("Variable site:", i+1)
                            #VariableSites += 1
                            Variable = True
                            break
                        #end if
                    #end if
                #end for
            #end with
            if Variable == False:
                InvariantSites += 1
            #end if
        #end for
      return InvariantSites
  #end method
  
  def Get_NTFrequencies(self, _Alignment, File_Type):
      count_A = 0
      count_T = 0
      count_C = 0
      count_G = 0
      count_N = 0
      total_num_nt = (self.NumCodonSites * 3) * self.NumSequences
      print("# Checking NT Frequencies")
      # Loop over sequences
      with open(_Alignment, "r") as handle:
          for record in tqdm(SeqIO.parse(handle, File_Type)):
              #print(record.id)
              for nt in str(record.seq):
                  if nt.upper() == "A": 
                      count_A += 1
                  if nt.upper() == "T": 
                      count_T += 1
                  if nt.upper() == "C": 
                      count_C += 1
                  if nt.upper() == "G": 
                      count_G += 1
                  if nt.upper() == "N": 
                      count_N += 1
                  #total_num_nt += 1
            #end for
          #end for
      #end with
      pi_A = count_A / total_num_nt
      pi_T = count_T / total_num_nt
      pi_C = count_C / total_num_nt
      pi_G = count_G / total_num_nt
      pi_N = count_N / total_num_nt
      return {"NucleotideFrequencies": {"A": pi_A,
                                        "T": pi_T,
                                        "G": pi_G,
                                        "C": pi_C,
                                        "N": pi_N
                                        }
              }
     
#end class

# =============================================================================
# Helper functions
# =============================================================================
"""
# =============================================================================
# Main
# =============================================================================
print()
print("#", "".join(["="]*78))
print("# Starting to profile multiple sequence alignment:", args.input)
print("#", "".join(["="]*78))
print()

if os.path.exists(args.input) and os.path.getsize(args.input) > 0:
    _MSA = AlignmentProfiler(args.input)
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
#https://www.promega.com/resources/guides/nucleic-acid-analysis/restriction-enzyme-resource/restriction-enzyme-resource-tables/iupac-ambiguity-codes-for-nucleotide-degeneracy/
print("# The alignment the following nucleotide frequencies...")
print()
print(tabulate([
    ['Adenine (A)',  _MSA.NTFrequencies["NucleotideFrequencies"]["A"]], 
    ['Thymine (T)',  _MSA.NTFrequencies["NucleotideFrequencies"]["T"]], 
    ['Guanine (G)',  _MSA.NTFrequencies["NucleotideFrequencies"]["G"]], 
    ['Cytosine (C)', _MSA.NTFrequencies["NucleotideFrequencies"]["C"]],
    ['Any Nucleotide (N)', _MSA.NTFrequencies["NucleotideFrequencies"]["N"]]], 
    headers=['Nucleotide', 'Frequency']))

#GC Content
#print()
#print("# The alignment the following GC content (%):", 
#      (_MSA.NTFrequencies["NucleotideFrequencies"]["G"] + 
#      _MSA.NTFrequencies["NucleotideFrequencies"]["C"]) * 100)

# Codon Frequencies
#print()
#print("# The alignment the following Codon frequencies")
#print()
#print(tabulate([
#    ['ATG',  0]
#    ], 
#    headers=['Codon', 'Frequency']))

# Search for ambiguous nucleotides
# Or ambiguous codons

# Pairwise Hamming distance
# Pairiwse p-distance
# Pairwise p-distance with indels
# Pairwise JC69 distance
# RSCU
# Pairwise KaKs

print()
"""
# =============================================================================
# Output CSV
# =============================================================================

"""
df = pd.DataFrame.from_dict(ResultsDict, orient='index')
print()
print("# Saving results to:", OutputCSV)
df.to_csv(OutputCSV, index=False)
"""



"""
# Summary stats
print("# Loaded an alignment from:", input_file)
print("# This alignment contains", TEST.num_sequences, "sequences")
print("# This alignment contains", TEST.num_sites, "sites")
print("# Number of invariant sites in alignment:", len(TEST.invariant_sites()))
print("# Fraction of invariant sites in alignment:", len(TEST.invariant_sites())/TEST.num_sites)
gap_list = TEST.gaps_distribution()
avg_gap = sum(gap_list)/ len(gap_list)
print("# Average measure of gappiness in alignment is:", avg_gap)
print("# Maximum measure of gappiness in alignment at a particular site is:", max(gap_list))
print("# Minimum measure of gappiness in alignment at a particular site is:", min(gap_list))

"""


# =============================================================================
# End of file
# =============================================================================
