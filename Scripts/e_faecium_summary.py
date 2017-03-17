#!/usr/bin/env python3

import itertools
import os
import copy
import csv
import re
import pyfastaq
import argparse

ariba_regex = re.compile(r'''(^VanB_Gly)|(^Van[HRSWXY]_B_Gly)''')
srst2_regex = re.compile(r'''(^VanB_Gly)|(^Van[HRSWXY]-B_Gly)''')
kmerres_regex = re.compile(r'''(^VanB)|(^Van[HRSWXY]-B)''')
read_depths = list(range(1,26,1)) + [30, 35, 40, 45, 50, 75, 100]


def load_suppl_data(infile):
    data = {}
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            assert 'SRR' in row
            data[row['SRR']] = row
    return fieldnames, data


def get_coverage_from_ariba_report(report_file, sample_dict):
    wanted_clusters = {x['ariba_cluster'] for x in sample_dict.values()}
    coverages = {x:[] for x in wanted_clusters}

    with open(report_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['cluster'] in wanted_clusters:
                coverages[row['cluster']].append(float(row['ctg_cov']))

    for x in coverages:
        if len(coverages[x]) == 0:
            coverages[x] = 0
        else:
            coverages[x] = round(sum(coverages[x]) / len(coverages[x]), 2)

    return coverages


def get_sequence_types(srr):
    with open(os.path.join(srr, 'ariba.mlst_report.tsv')) as f:
        f.readline()
        ariba_st = f.readline().split()[0]

    with open(os.path.join(srr, 'srst2.mlst_results.txt')) as f:
        f.readline()
        srst2_st = f.readline().split()[1]

    return ariba_st, srst2_st


def load_ariba_csv(infile):
    header = None
    data = {}

    for line in open(infile):
        if header is None:
            header = line.rstrip().split(',')
        else:
            fields = line.rstrip().split(',')
            assert len(header) == len(fields)
            d = {x[0]: x[1] for x in zip(header[1:], fields[1:]) if ariba_regex.match(x[0])}
            new_d = {}
            for key, val in d.items():
                group = key.split('_')[0]
                ariba_cluster, match_or_ref = key.split('.')
                if group not in new_d:
                    new_d[group] = {'match': None, 'ref_seq': None, 'ariba_cluster': ariba_cluster}
                if match_or_ref == 'ref_seq':
                    if val != 'NA':
                        cluster = val.split('.')[0][:-4].replace('_', '-')
                        ref = val.split('_')[-1]
                        val = cluster + '_' + ref
                new_d[group][match_or_ref] = val

            data[fields[0]] = new_d

    return data


def load_ariba_results(rootdir, outprefix, samples, data, subsample_data):
    ariba_summary_prefix = outprefix + '.ariba.summary'
    pyfastaq.utils.syscall('ariba summary --no_tree --cluster_cols assembled,ref_seq ' + ariba_summary_prefix + ' ' + rootdir + '/*/ariba.report.tsv ' + rootdir + '/*/Sample_reads/*/ariba.report.tsv')
    raw_data = load_ariba_csv(ariba_summary_prefix + '.csv')

    for sample in samples:
        sample_filename = os.path.join(rootdir, sample, 'ariba.report.tsv')
        assert sample_filename in raw_data
        assert sample in data
        data[sample]['ARIBA'] = raw_data[sample_filename]
        coverages = get_coverage_from_ariba_report(sample_filename, data[sample]['ARIBA'])
        for d in data[sample]['ARIBA'].values():
            d['mean_coverage'] = coverages[d['ariba_cluster']]

        if samples[sample]:
            assert sample in subsample_data
            for i in read_depths:
                sample_filename = os.path.join(rootdir, sample, 'Sample_reads', str(i), 'ariba.report.tsv')
                assert sample_filename in raw_data
                subsample_data[sample][i]['ARIBA'] = raw_data[sample_filename]



def load_srst2_file(infile):
    with open(infile) as f:
        header = f.readline().rstrip().split()
        fields = f.readline().rstrip().split()
        assert len(header) == len(fields)

    data = {}

    for key, val in zip(header, fields):
        if srst2_regex.match(key):
            key = key.split('_')[0].split('-')[0]
            d = {'?': False, '*': False}
            for c in ['?', '*']:
                if val.endswith(c):
                    d[c] = True
                    val = val.rstrip(c)

            d['ref_seq'] = val
            data[key] = d

    return data


def add_depths_from_fullgenes_file(data_dict, infile):
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames

        for row in reader:
            short_gene = row['gene'].split('_')[0].split('-')[0]
            if short_gene in data_dict:
                assert data_dict[short_gene]['ref_seq'] == row['allele']
                data_dict[short_gene]['depth'] = row['depth']


def load_srst2_results(rootdir, samples, data, subsample_data):
    for sample in samples:
        sample_filename = os.path.join(rootdir, sample, 'srst2.out__genes__results.txt')
        sample_fullgenes_filename = os.path.join(rootdir, sample, 'srst2.out__fullgenes__results.txt')
        data[sample]['SRST2'] = load_srst2_file(sample_filename)
        add_depths_from_fullgenes_file(data[sample]['SRST2'], sample_fullgenes_filename)

        if samples[sample]:
            assert sample in subsample_data

            for i in read_depths:
                sample_filename = os.path.join(rootdir, sample, 'Sample_reads', str(i), 'srst2.out__genes__results.txt')
                subsample_data[sample][i]['SRST2'] = load_srst2_file(sample_filename)


def load_kmerres_file(infile):
    hits = {}

    for line in open(infile):
        data = line.rstrip().split()
        final_col_split = data[-1].split(';')
        if ';' in data[-1] and kmerres_regex.match(final_col_split[2]) and final_col_split[3] == 'Gly':
            cluster = final_col_split[2].split('-')[0]
            ref_seq = final_col_split[2] + '_' + data[1]
            hits[cluster] = {'ref_seq': ref_seq, 'pval': data[10], 'coverage': data[5], 'depth': data[7]}

    return hits


def load_kmerres_results(rootdir, samples, data, subsample_data):
    for sample in samples:
        sample_filename = os.path.join(rootdir, sample, 'kmerResistance.out1')
        data[sample]['KmerResistance'] = load_kmerres_file(sample_filename)

        if samples[sample]:
            assert sample in subsample_data

            for i in read_depths:
                sample_filename = os.path.join(rootdir, sample, 'Sample_reads', str(i), 'kmerResistance.out1')
                subsample_data[sample][i]['KmerResistance'] = load_kmerres_file(sample_filename)



parser = argparse.ArgumentParser(
    description = 'Gathers vanB data',
)

parser.add_argument('rootdir', help='E_faecium/ directory')
parser.add_argument('outprefix', help='Prefix of output files')
parser.add_argument('--lenient', action='store_true', help='Also count SRST2 with ? as present, and ariba partial/fragmented/interrupted as present')
options = parser.parse_args()

print('load supplementary data')
suppl_data_columns, suppl_data = load_suppl_data('howden_supplementary_with_srr.tsv')

samples = [x for x in os.listdir(options.rootdir) if x.startswith('SRR')]
sample0 = samples[0]
samples = {x: os.path.exists(os.path.join(options.rootdir, x, 'Sample_reads')) for x in samples}
data = {x: {} for x in samples}
subsample_data = {x: {i: {} for i in read_depths} for x in samples if samples[x]}

print('load ARIBA results', flush=True)
load_ariba_results(options.rootdir, options.outprefix, samples, data, subsample_data)
print('load SRST2 results', flush=True)
load_srst2_results(options.rootdir, samples, data, subsample_data)
print('load KmerResistance results', flush=True)
load_kmerres_results(options.rootdir, samples, data, subsample_data)

f = open(options.outprefix + '.supplementary_table.tsv', 'w')
sample_list = sorted(list(samples.keys()))
tools = ['ARIBA', 'KmerResistance', 'SRST2']
clusters = sorted(list(data[sample0]['ARIBA'].keys()))

column_headers = copy.copy(suppl_data_columns)
column_headers.extend(['ARIBA_ST', 'SRST2_ST'])
for cluster in clusters:
    for tool in tools:
        for suffix in ['', '.ref_seq', '.depth']:
            column_headers.append(tool + '.' + cluster + suffix)

print(*column_headers, sep='\t', file=f)

for sample in sorted(data):
    line = [suppl_data[sample][x] for x in suppl_data_columns]
    ariba_st, srst2_st = get_sequence_types(sample)
    line.extend([ariba_st, srst2_st])

    for cluster in sorted(clusters):
        line.append(data[sample]['ARIBA'][cluster]['assembled'])
        line.append(data[sample]['ARIBA'][cluster]['ref_seq'])
        line.append(data[sample]['ARIBA'][cluster]['mean_coverage'])

        if cluster in data[sample]['KmerResistance']:
            line.append(data[sample]['KmerResistance'][cluster]['coverage'])
            line.append(data[sample]['KmerResistance'][cluster]['ref_seq'])
            line.append(data[sample]['KmerResistance'][cluster]['depth'])
        else:
            line.append('0')
            line.append('NA')
            line.append('NA')

        if cluster in data[sample]['SRST2']:
            d = data[sample]['SRST2'][cluster]
            line.append('yes')
            if d['*']:
                line[-1] += '*'
            if d['?']:
                line[-1] += '?'
            line.append(d['ref_seq'])
            line.append(d['depth'])
        else:
            line.append('no')
            line.append('NA')
            line.append('NA')

    print(*line, sep='\t', file=f)


f.close()

skip_clusters = {'VanW', 'VanY'}
total_calls_per_depth = {x: {'ARIBA': 0, 'SRST2': 0, 'KmerResistance': 0} for x in read_depths}
r_script = options.outprefix + '.make_plots.R'
r_fh = pyfastaq.utils.open_file_write(r_script)
print('library(ggplot2)', file=r_fh)
print('library(grid)', file=r_fh)
print(r'''plot_style = theme_bw()+
    theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(size=16),
    axis.title.y = element_text(size=20),
    axis.text.y = element_text(size=16),
    axis.title = element_text(size=20),
    legend.title = element_text(size=18),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(lineheight=.6, size = 24))

plot_style_no_legend = theme_bw()+
    theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    #axis.title.x = element_blank(),
    axis.title.y = element_text(size=20),
    #axis.text.x = element_blank(),
    axis.text.y = element_text(size=20),
    axis.title = element_text(size=20),
    plot.title = element_text(lineheight=.6, size = 24),
    legend.position="none")
''', file=r_fh)

for cluster in clusters:
    if cluster in skip_clusters:
        continue

    print('Make plot data for cluster', cluster)
    outprefix = options.outprefix + '.depth_plot.' + cluster
    csv_file = outprefix + '.csv'
    with open(csv_file, 'w') as csv_fh:
        print('Read_depth,Tool,Number_of_calls', file=csv_fh)

        for i in read_depths:
            calls_per_depth = {'ARIBA': 0, 'SRST2': 0, 'KmerResistance': 0}
            for sample in subsample_data:
                ariba_assembled = subsample_data[sample][i]['ARIBA'][cluster]['assembled']

                if ariba_assembled.startswith('yes') or (options.lenient and ariba_assembled != 'no'):
                    calls_per_depth['ARIBA'] += 1

                d = subsample_data[sample][i]['SRST2']
                if cluster in d and (options.lenient or not d[cluster]['?']):
                    calls_per_depth['SRST2'] += 1

                d = subsample_data[sample][i]['KmerResistance']
                if cluster in d:
                    calls_per_depth['KmerResistance'] += 1


            for tool in calls_per_depth:
                print(i, tool, calls_per_depth[tool], sep=',', file=csv_fh)
                total_calls_per_depth[i][tool] += calls_per_depth[tool]

    print('d=read.csv(file="', csv_file, '", header=TRUE)', sep='', file=r_fh)
    print('ggplot(data=d) +', file=r_fh)
    print('    geom_hline(yintercept=', len(subsample_data), ', linetype=2) +', sep='', file=r_fh)
    print('    geom_line(mapping=aes(x=Read_depth, y=Number_of_calls, group=Tool, color=Tool)) +', file=r_fh)
    print('    geom_point(mapping=aes(x=Read_depth, y=Number_of_calls, group=Tool, color=Tool, shape=Tool)) +', file=r_fh)
    print('    xlab("Mean read depth") +', file=r_fh)
    print('    ylab("Number of calls") +', file=r_fh)
    print('    plot_style +', file=r_fh)
    print('    scale_color_brewer(palette = "Accent") +', file=r_fh)
    print('ggsave(filename="', outprefix + '.pdf", width=8, height=3)', sep='', file=r_fh)

print('Make summary plot data')
outprefix = options.outprefix + '.depth_plot.all'
csv_file = outprefix + '.csv'
with open(csv_file, 'w') as f:
    print('Read_depth,Tool,Number_of_calls', file=f)
    for i in total_calls_per_depth:
        for tool in total_calls_per_depth[i]:
            print(i, tool, total_calls_per_depth[i][tool], sep=',', file=f)

print('d=read.csv(file="', csv_file, '", header=TRUE)', sep='', file=r_fh)
print('ggplot(data=d) +', file=r_fh)
print('    geom_hline(yintercept=', len(subsample_data) * (len(clusters) - len(skip_clusters)), ', linetype=2) +', sep='', file=r_fh)
print('    geom_line(mapping=aes(x=Read_depth, y=Number_of_calls, group=Tool, color=Tool)) +', file=r_fh)
print('    geom_point(mapping=aes(x=Read_depth, y=Number_of_calls, group=Tool, color=Tool, shape=Tool)) +', file=r_fh)
print('    xlab("Mean read depth") +', file=r_fh)
print('    ylab("Number of calls") +', file=r_fh)
print('    plot_style +', file=r_fh)
print('    scale_color_brewer(palette = "Accent") +', file=r_fh)
print('ggsave(filename="', outprefix + '.pdf", width=8, height=3)', sep='', file=r_fh)

pyfastaq.utils.close(r_fh)
print('Run R script to make plots')
pyfastaq.utils.syscall('R CMD BATCH ' + r_script)
try:
    os.unlink('Rplots.pdf')
except:
    pass
