import csv
import sys
import logging

"""
Given a set of Marcin's biclusters, map known and clinically validated syn lethal
gene pairs to biclusters
"""

ncbi_to_symbol_file = "data/human.gene_info"

sym_to_ncbi = {}
ncbi_to_sym = {}

# make ncbi id to gene symbol map
with open(ncbi_to_symbol_file) as ns_file:
    reader = csv.reader(ns_file, delimiter="\t")
    for row in reader:
        sym_to_ncbi[row[2]] = row[1]
        ncbi_to_sym[row[1]] = row[2]

sys.argv.pop(0)
for gene in sys.argv:
    gene = gene.strip()
    if gene in sym_to_ncbi:
        print(gene + ": NCBI:" + sym_to_ncbi[gene])
    else:
        logging.warning("Couldn't find ncbi id for " + gene)

