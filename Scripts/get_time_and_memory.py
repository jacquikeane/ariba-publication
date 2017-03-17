#!/usr/bin/env python3
import os
import sys

def read_bash_out_file(infile):
    data = {}
    with open(infile) as f:
        for line in f:
            if 'wall clock' in line:
                fields = line.rstrip().split()
                time_list = fields[-1].split(':')
                if len(time_list) == 2:
                    seconds = 60 * int(time_list[0]) + int(round(float(time_list[1]), 0))
                else:
                    seconds = 60 * 60 * int(time_list[0]) + 60 * int(time_list[1]) + int(time_list[2])
                if fields[0] not in data:
                    data[fields[0]] = {}
                data[fields[0]]['time'] = seconds
            elif 'Maximum resident set size' in line:
                fields = line.rstrip().split()
                kb = int(fields[-1])
                if fields[0] not in data:
                    data[fields[0]] = {}
                data[fields[0]]['memory'] = round(kb / 1024, 2)

    return data


def write_resources_file(data, outfile):
    with open(outfile, 'w') as f:
        print('Dataset\tSample\tTool\tWall_clock\tRAM', file=f)

        for dataset in data:
            for sample in data[dataset]['ARIBA']:
                for tool in data[dataset]:
                    print(dataset, sample, tool, data[dataset][tool][sample]['time'], data[dataset][tool][sample]['memory'], sep='\t', file=f)


root_dir = sys.argv[1]


data = {
    'E.faecium.mlst': {
        'ARIBA': read_bash_out_file(os.path.join(root_dir, 'E_faecium', 'bash_out.ariba.mlst')),
        'SRST2': read_bash_out_file(os.path.join(root_dir, 'E_faecium', 'bash_out.srst2.mlst')),
    },
    'E.faecium': {
        'ARIBA': read_bash_out_file(os.path.join(root_dir, 'E_faecium', 'bash_out.ariba')),
        'KmerResistance': read_bash_out_file(os.path.join(root_dir, 'E_faecium', 'bash_out.kmerResistance')),
        'SRST2': read_bash_out_file(os.path.join(root_dir, 'E_faecium', 'bash_out.srst2')),
    },
    'S.sonnei': {
        'ARIBA': read_bash_out_file(os.path.join(root_dir, 'S_sonnei', 'bash_out.ariba')),
        'KmerResistance': read_bash_out_file(os.path.join(root_dir, 'S_sonnei', 'bash_out.kmerResistance')),
        'SRST2': read_bash_out_file(os.path.join(root_dir, 'S_sonnei', 'bash_out.srst2')),
    }
}


write_resources_file(data, os.path.join(root_dir, 'Time_and_memory', 'resources.tsv'))

