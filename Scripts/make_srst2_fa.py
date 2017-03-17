#!/usr/bin/env python3

import os
import argparse
import pyfastaq

parser = argparse.ArgumentParser(
    description = 'Makes SRST2 input file from ariba db directory',
    usage = '%(prog)s <prepareref dir> <outprefix>'
)

parser.add_argument('ariba_db_dir', help='ARIBA prepareref dir')
parser.add_argument('outfile', help='Name of output fasta')
options = parser.parse_args()

ariba_all_seqs_fa = os.path.join(options.ariba_db_dir, '02.cdhit.all.fa')
ariba_clusters_file = os.path.join(options.ariba_db_dir, '02.cdhit.clusters.tsv')

seq2clusters = {}
cluster2id = {}

for line in open(ariba_clusters_file):
    cluster_name, *seq_names = line.rstrip().split()
    if cluster_name not in cluster2id:
        cluster2id[cluster_name] = len(cluster2id) + 1

    for seq_name in seq_names:
        assert seq_name not in seq2clusters
        seq2clusters[seq_name] = cluster_name


seq_counter = 1
f = pyfastaq.utils.open_file_write(options.outfile)

for seq in pyfastaq.sequences.file_reader(ariba_all_seqs_fa):
    cluster_name = seq2clusters[seq.id]
    cluster_id = cluster2id[cluster_name]
    srst2_name = '__'.join([
        str(cluster_id),
        cluster_name.replace('_', '-'),
        seq.id.replace('_', '-'),
        str(seq_counter)
    ])

    seq.id = srst2_name + ' ' + seq.id
    print(seq, file=f)
    seq_counter += 1

pyfastaq.utils.close(f)

