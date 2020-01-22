import csv
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter

"""
Given a set of Marcin's biclusters, map known and clinically validated syn lethal
gene pairs to biclusters
"""

bicluster_file = "data/results_gene_effect_corrected__nr_0.25_score_root_transpose.txt"
x_coordinate_to_ncbi_file = "data/gene_effect_corrected_transpose_NCBI_geneids.txt"
ncbi_to_symbol_file = "data/human.gene_info"
syn_lethal_file = "data/gene_gene_sl_relationships.csv"

sl_pairs = []
sym_to_ncbi = {}
ncbi_to_sym = {}

x_coordinate_to_ncbi = {}

bicluster_to_ncbi = {}
ncbi_to_bicluster = {}

# read in syn lethal gene pairs
with open(syn_lethal_file) as csv_file:
    reader = csv.DictReader(filter(lambda row: row[0] != '#', csv_file))
    line_count = 0
    for row in reader:
        sl_pairs.append(row)
        line_count += 1

# make ncbi id to gene symbol map
with open(ncbi_to_symbol_file) as ns_file:
    reader = csv.reader(ns_file, delimiter="\t")
    for row in reader:
        sym_to_ncbi[row[2]] = row[1]
        ncbi_to_sym[row[1]] = row[2]

# make of x coordinate to ncbi gene
with open(x_coordinate_to_ncbi_file) as cn_file:
    reader = csv.reader(cn_file, delimiter="\t")
    for row in reader:
        x_coordinate_to_ncbi[row[0]] = row[1]

# make cluster to gene links
with open(bicluster_file) as bc_file:
    reader = csv.DictReader(bc_file, delimiter="\t")
    for row in reader:
        for x_coord in row['block_id'].split("/")[0].split(","):
            ncbi_id = x_coordinate_to_ncbi[x_coord]

            if not row['#number'] in bicluster_to_ncbi:
                bicluster_to_ncbi[row['#number']] = []
            bicluster_to_ncbi[row['#number']].append(ncbi_id)

            if not ncbi_id in ncbi_to_bicluster:
                ncbi_to_bicluster[ncbi_id] = []
            ncbi_to_bicluster[ncbi_id].append(row['#number'])

edges = []
G = nx.Graph()

# make edges
for dict_entry in sl_pairs:
    if dict_entry['SL'] == '0':
        continue
    try:
        gene1_biclusters = ncbi_to_bicluster["NCBI:" + sym_to_ncbi[dict_entry['gene1']]]
    except KeyError as ke:
        print("problem getting gene 1 biclusters: " + str(ke))
    try:
        gene2_biclusters = ncbi_to_bicluster["NCBI:" + sym_to_ncbi[dict_entry['gene2']]]
    except KeyError as ke:
        print("problem getting gene 2 biclusters: " + str(ke))
    for b1 in gene1_biclusters:
        for b2 in gene2_biclusters:
            G.add_node(b1)
            G.add_node(b2)
            if G.has_edge(b1, b2):
                # we added this one before, just increase the weight by one
                G[b1][b2]['weight'] += 1
            else:
                # new edge. add with weight=1
                G.add_edge(b1, b2, weight=1)
            # G.add_edge(b1, b2)
            edges.append([b1, b2])

plt.subplot(121)
nx.draw(G, with_labels=False)
plt.show()
pr = nx.pagerank(G)
sorted_nodes = {k: v for k, v in sorted(pr.items(), key=lambda item: item[1])}

w = csv.writer(open("nodes_by_page_rank.csv", "w"))
for key, val in sorted_nodes.items():
    w.writerow([key, val])

# write self connected nodes, by weight
self_conn_nodes = []
for edge in G.edges(data=True):
    if edge[0] == edge[1]:
        self_conn_nodes.append(edge)

self_conn_nodes = sorted(
    self_conn_nodes, key=lambda item: item[2]['weight'], reverse=True)
w = csv.writer(open("self_connected_nodes.csv", "w"))
for node in self_conn_nodes:
    w.writerow(node)