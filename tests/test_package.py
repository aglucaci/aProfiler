#Primate_TP53.codons.cln.fa

clear

testData="data/Primate_TP53.codons.cln.fa"

echo "[INFO] Running tests on: $testData"


# Test basic analysis and entropy plot
alignmentprofiler data/Primate_TP53.codons.cln.fa --format fasta --output summary.json --plot entropy --plot-out entropy.png

# Test gaps plot
alignmentprofiler data/Primate_TP53.codons.cln.fa --format fasta --output summary.json --plot gaps --plot-out gaps.png

# Test codon plot
alignmentprofiler data/Primate_TP53.codons.cln.fa --format fasta --output summary.json --plot codon_diversity --plot-out codon_diversity.png

# Test codon plot
alignmentprofiler data/Primate_TP53.codons.cln.fa --format fasta --output summary.json --plot protein_diversity --plot-out protein_diversity.png


# Test csv output
alignmentprofiler data/Primate_TP53.codons.cln.fa --output summary.json --csv sequences.csv

# Test RSCU file
alignmentprofiler $testData --rscu-csv rscu_table.csv

# Test RSCU file
alignmentprofiler $testData --aa-csv aa_table.csv


exit 0
