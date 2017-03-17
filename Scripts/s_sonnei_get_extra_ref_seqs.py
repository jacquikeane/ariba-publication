#!/usr/bin/env python3

import os
import sys
import shutil
import urllib.request
import tempfile
import pyfastaq

def load_ref_data(infile):
    data = {}

    with open(infile) as f:
        for line in f:
            ref_name, protein_id, start, end, antibio = line.rstrip().split('\t')
            if ref_name == 'Reference_ID':
                continue
            if ref_name not in data:
                data[ref_name] = set()
            data[ref_name].add((protein_id, int(start) - 1, int(end) - 1, antibio))

    return data


def get_ref_seq(ref_name):
    tmpdir = tempfile.mkdtemp(prefix='tmp.shigella_get_extra_ref_seqs.' + ref_name + '.', dir=os.getcwd())
    fasta_file = os.path.join(tmpdir, 'seq.fa')
    print('Getting', ref_name, flush=True)
    try:
        urllib.request.urlretrieve('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=' + ref_name, filename=fasta_file)
    except:
        shutil.rmtree(tmpdir)
        print('Error getting sequence', ref_name, file=sys.stderr)

    print('   ... got', ref_name, flush=True)
    seqs = {}
    pyfastaq.tasks.file_to_dict(fasta_file, seqs)
    assert len(seqs) == 1
    seq_name, ref_seq = seqs.popitem()
    shutil.rmtree(tmpdir)
    return ref_seq


ref_data = load_ref_data('ref_data.in.tsv')
f_fa = pyfastaq.utils.open_file_write('ref_data.sequences.fa')
f_tsv = pyfastaq.utils.open_file_write('ref_data.metadata.tsv')

for ref_name, protein_set in ref_data.items():
    ref_seq = get_ref_seq(ref_name)

    for (protein_name, start, end, antibio) in protein_set:
        new_seq = ref_seq.subseq(start, end+1)
        original_name = new_seq.id
        new_seq.id = antibio + '.' + '_'.join(['manually_added', ref_name, protein_name, str(start + 1), str(end + 1)])
        print(new_seq.id, '1', '0', '.', '.', 'Original ref seq: ' + original_name, sep='\t', file=f_tsv)
        print(new_seq, file=f_fa)
        print('Made sequence', new_seq.id)

pyfastaq.utils.close(f_fa)
pyfastaq.utils.close(f_tsv)

