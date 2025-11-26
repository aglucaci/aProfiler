# alignmentprofiler/io.py
from Bio import AlignIO

def load_alignment(filepath, format='fasta'):
    return AlignIO.read(filepath, format)
# end method


