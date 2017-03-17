#!/usr/bin/env python3

import math
import pprint
import copy
import csv
import sys
import os
import operator
import argparse
import pyfastaq


tool_colours = {
    'ariba': '#7fc97f',
    'kmerResistance': '#beaed4',
    'srst2': '#fdc086',
    'Ref': '#ffff99',
}


class Sample:
    def __init__(self, ref_data, lenient):
        self.ref_data = ref_data
        self.name = ref_data['ERR']
        self.min_ref_genotype_cov = 90
        self.ref_genotypes = Sample._get_genotypes(ref_data, self.min_ref_genotype_cov)
        self.ariba_raw_data = None
        self.tool_data = {'ariba': {}, 'srst2': {}, 'kmerResistance': {}}
        self.antibio_calls = {tool: {k: None for k in self.ref_genotypes} for tool in self.tool_data}
        self.antibio_depths = {tool: {k: None for k in self.ref_genotypes} for tool in self.tool_data}
        self.ariba_report_data = {} # cluster => list of report dicts
        self.lenient = lenient


    def __str__(self):
        lines = []
        for key in self.ref_genotypes:
            lines.append('\t'.join([
                self.name,
                'calls(geno,pheno,ariba,kmerres,srst2)',
                key,
                str(self.ref_genotypes[key]),
                self.ref_data[key],
                str(self.antibio_calls['ariba'][key]),
                str(self.antibio_calls['kmerResistance'][key]),
                str(self.antibio_calls['srst2'][key])
            ]))

        for tool in ['ariba', 'kmerResistance', 'srst2']:
            for key in sorted(self.tool_data[tool]):
                lines.append('\t'.join([self.name , tool, key, str(self.tool_data[tool][key])]))

        return '\n'.join(lines)


    @classmethod
    def _get_one_genotype(cls, ref_data, genes, min_cov):
        for gene in genes:
            if float(ref_data[gene]) >= min_cov:
                return True
        return False


    @classmethod
    def _get_genotypes(cls, ref_data, min_cov):
        genotypes = {}

        has_strA = Sample._get_one_genotype(ref_data, ['strA'], min_cov)
        has_strB = Sample._get_one_genotype(ref_data, ['strB'], min_cov)
        has_str = Sample._get_one_genotype(ref_data, ['aadA1', 'aadA2'], min_cov) or (has_strA and has_strB)

        has_tet1 = Sample._get_one_genotype(ref_data, ['tetA', 'tetR'], min_cov)
        has_tet2 = Sample._get_one_genotype(ref_data, ['tetA.1', 'tetC', 'tetD', 'tetR.1'], min_cov)

        genotypes = {
            'Amp': Sample._get_one_genotype(ref_data, ['blaTEM', 'bla.OXA.1', 'bla.OXA.10', 'bla.CTXM.15'], min_cov),
            'Cmp': Sample._get_one_genotype(ref_data, ['catA1', 'catB3'], min_cov),
            'GyrA': ref_data['GyrA'] != 'no',
            'Str': has_str,
            'Sul': Sample._get_one_genotype(ref_data, ['sul1', 'sul2'], min_cov),
            'Tet': has_tet1 or has_tet2,
            'Tri': Sample._get_one_genotype(ref_data, ['dfrA1', 'dhfr3b', 'dfrA5', 'dfrA8', 'dfrA14', 'dfrA12'], min_cov),
        }

        return genotypes


    def load_kmerres_file(self, infile):
        for line in open(infile):
            if line.startswith('#Template'):
                continue
            data = line.rstrip().split()
            self.tool_data['kmerResistance'][data[0].replace('-', '_')] = float(data[7])


    @classmethod
    def _load_srst2_fullgenes_file(cls, infile):
        data = {}
        with open(infile) as f:
            for line in f:
                if line.startswith('Sample'):
                    continue

                fields = line.rstrip().split('\t')
                data[fields[3]] = fields[-1], float(fields[5])

        return data


    def load_srst2_genes_file(self, fullgenes_file, genes_file):
        name_dict = Sample._load_srst2_fullgenes_file(fullgenes_file)

        with open(genes_file) as f:
            header = f.readline().rstrip().split()
            fields = f.readline().rstrip().split()
            assert len(header) == len(fields)

        for key, val in zip(header[1:], fields[1:]):
            d = {'?': False, '*': False}
            for c in ['?', '*']:
                if val.endswith(c):
                    d[c] = True
                    val = val.rstrip(c)

            self.tool_data['srst2'][name_dict[val][0]] = d
            d['depth'] = name_dict[val][1]


    def load_ariba_report(self, report_tsv):
        with open(report_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['cluster'] not in self.ariba_report_data:
                    self.ariba_report_data[row['cluster']] = []
                self.ariba_report_data[row['cluster']].append(copy.copy(row))


    def _has_gene(self, tool, gene, cluster=None):
        depth = None
        if tool == 'ariba':
            assert cluster is not None
            assembled = self.tool_data['ariba'].get(cluster + '.assembled', 'no')
            assembled_ok = assembled.startswith('yes') or (self.lenient and assembled != 'no')
            has_gene = assembled_ok and self.tool_data['ariba'].get(cluster + '.ref_seq', None) == gene

            if has_gene:
                coverages = [float(x['ctg_cov']) for x in self.ariba_report_data[cluster]]
                if len(coverages):
                    depth = round(sum(coverages) / len(coverages), 1)
        elif tool == 'srst2':
            has_gene = gene in self.tool_data['srst2'] and (self.lenient or not self.tool_data['srst2'][gene]['?'])
            if has_gene:
                depth = self.tool_data['srst2'][gene]['depth']
        elif tool == 'kmerResistance':
            has_gene = gene in self.tool_data['kmerResistance']
            if has_gene:
                depth = self.tool_data['kmerResistance'][gene]

        return has_gene, depth


    def _has_any_gene_from_set_of_clusters(self, tool, gene_set, cluster=None):
        if tool == 'ariba':
            assert cluster is not None

        for gene in gene_set:
            has_gene, depth = self._has_gene(tool, gene, cluster=cluster)
            if has_gene:
                return True, depth
        return False, None


    def _has_any_cluster_of_genes(self, tool, genes_dict):
        for cluster, gene_set in genes_dict.items():
            has_any_gene, depth = self._has_any_gene_from_set_of_clusters(tool, gene_set, cluster=cluster)
            if has_any_gene:
                return True, depth
        return False, None


    def _call_gyra(self):
        for tool in self.antibio_calls:
            self.antibio_calls[tool]['GyrA'] = False

        cluster = 'gyrA_7'

        if self._has_gene('ariba', 'gyrA.3003294.U00096.2336792_2339420.2052', cluster=cluster):
            for snp in ['S83L', 'D87Y', 'D87G']:
                if self.tool_data['ariba'].get(cluster + '.' + snp, None) == 'yes':
                    self.antibio_calls['ariba']['GyrA'] = True
                    coverages = [float(x['ctg_cov']) for x in self.ariba_report_data[cluster]]
                    if len(coverages):
                        self.antibio_depths['ariba']['GyrA'] = round(sum(coverages) / len(coverages), 1)
                    break


    def _call_amp_sul_tri_cmp(self, resistance_clusters):
        for antibio in {'Amp', 'Sul', 'Tri', 'Cmp'}:
            for tool in self.antibio_calls:
                self.antibio_calls[tool][antibio], self.antibio_depths[tool][antibio] = self._has_any_cluster_of_genes(tool, resistance_clusters[antibio])



    def _call_str(self, resistance_clusters):
        for tool in self.antibio_calls:
            self.antibio_calls[tool]['Str'] = False

            has_Str, Str_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Str'])
            if has_Str:
                self.antibio_calls[tool]['Str'] = True
                self.antibio_depths[tool]['Str'] = Str_depth
                continue

            has_Str_A, Str_A_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Str_A'])
            has_Str_B, Str_B_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Str_B'])
            if has_Str_A and has_Str_B:
                self.antibio_calls[tool]['Str'] = True
                if None not in [has_Str_A, has_Str_B]:
                    self.antibio_depths[tool]['Str'] = 0.5 * (Str_A_depth + Str_B_depth)


    def _call_tet(self, resistance_clusters):
        for tool in self.antibio_calls:
            self.antibio_calls[tool]['Tet'] = False

            has_Tet_AR_A, Tet_AR_A_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Tet_AR_A'])
            has_Tet_AR_R, Tet_AR_R_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Tet_AR_R'])
            if has_Tet_AR_A and has_Tet_AR_R:
                self.antibio_calls[tool]['Tet'] = True
                if None not in [Tet_AR_A_depth, Tet_AR_R_depth]:
                    self.antibio_depths[tool]['Tet'] = 0.5 * (Tet_AR_A_depth + Tet_AR_R_depth)
                continue

            has_Tet_ACDR_A, Tet_ACDR_A_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Tet_ACDR_A'])
            has_Tet_ACDR_C, Tet_ACDR_C_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Tet_ACDR_C'])
            has_Tet_ACDR_D, Tet_ACDR_D_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Tet_ACDR_D'])
            has_Tet_ACDR_R, Tet_ACDR_R_depth = self._has_any_cluster_of_genes(tool, resistance_clusters['Tet_ACDR_R'])
            if has_Tet_ACDR_A and has_Tet_ACDR_C and has_Tet_ACDR_D and has_Tet_ACDR_R:
                self.antibio_calls[tool]['Tet'] = True
                if None not in [Tet_ACDR_A_depth, Tet_ACDR_C_depth, Tet_ACDR_D_depth, Tet_ACDR_R_depth]:
                    self.antibio_depths[tool]['Tet'] = 0.25 * (Tet_ACDR_A_depth + Tet_ACDR_C_depth + Tet_ACDR_D_depth + Tet_ACDR_R_depth)


    def call_resistance(self, resistance_clusters):
        self._call_gyra()
        self._call_amp_sul_tri_cmp(resistance_clusters)
        self._call_str(resistance_clusters)
        self._call_tet(resistance_clusters)


class Samples:
    def __init__(self, ref_data_tsv, srst2_ref_fasta, ariba_clusters_tsv, root_indir, files_outprefix, lenient):
        self.ref_data_tsv = ref_data_tsv
        self.srst2_ref_fasta = srst2_ref_fasta
        self.ariba_clusters_tsv = ariba_clusters_tsv
        self.rootdir = root_indir
        self.outprefix = files_outprefix
        self.samples = {}  # sample name -> Sample object
        self.tools = ['ariba', 'kmerResistance', 'srst2', 'Ref']
        self.lenient = lenient


    def _load_resistance_genes(self):
        resistance_clusters = {
            'Amp': {'TEM_1+', 'OXA_1+', 'OXA_10+', 'CTX_M_1+'},
            'Cmp': {'cat+', 'catB-'},
            'Str': {'aad-', 'aadA-'},
            'Str_A': {'APH_3____Ib+'},
            'Str_B': {'APH_6__Id+'},
            'Sul': {'sul1', 'sul2'},
            'Tet_AR_A': {'tetA'},
            'Tet_AR_R': {'tetR'},
            'Tet_ACDR_A': {'tet-'},
            'Tet_ACDR_C': {'tetC'},
            'Tet_ACDR_D': {'tetD'},
            'Tet_ACDR_R': {'tetR_1'},
            'Tri': {'dfrA1_1', 'dfrA1', 'dfhr3b', 'dfrA5', 'dfrA8', 'dfrA12', 'dfrA14'},
        }

        clusters = {}
        with open(self.ariba_clusters_tsv) as f:
            for line in f:
                data = line.rstrip().split('\t')
                assert data[0] not in clusters
                clusters[data[0]] = set(data[1:])

        self.resistance_genes = {}
        for antibio, cluster_set in resistance_clusters.items():
            self.resistance_genes[antibio] = {x: clusters[x] for x in cluster_set}


    def _write_amr_genes_suppl_tsv(self):
        ariba_to_srst2 = {}
        with open(self.srst2_ref_fasta) as f:
            for line in f:
                if line.startswith('>'):
                    srst2_name, ariba_name = line.rstrip()[1:].split()
                    assert ariba_name not in ariba_to_srst2
                    ariba_to_srst2[ariba_name] = srst2_name

        with open(self.outprefix + '.amr_genes.tsv', 'w') as f:
            print('Antibiotic', 'ARIBA cluster', 'ARIBA name', 'SRST2 name', sep='\t', file=f)
            for antibio, cluster_dict in sorted(self.resistance_genes.items()):
                for cluster, seq_set in sorted(list(cluster_dict.items())):
                    for ariba_name in sorted(seq_set):
                        srst2_name = ariba_to_srst2[ariba_name]
                        print(antibio[:3], cluster, ariba_name, srst2_name, sep='\t', file=f)

            gyra_seq = 'gyrA.3003294.U00096.2336792_2339420.2052'
            print('Nal', 'gyrA_7', gyra_seq, ariba_to_srst2[gyra_seq], sep='\t', file=f)


    def _load_ref_tsv(self):
        with open(self.ref_data_tsv) as f:
            header_dict = None

            for line in f:
                line_data = line.rstrip().split('\t')
                if header_dict is None:
                    header_dict = {line_data[i]: i for i in range(len(line_data))}
                else:
                    data_dict = {x: line_data[header_dict[x]] for x in header_dict}
                    sample = Sample(data_dict, self.lenient)
                    assert sample.name not in self.samples
                    self.samples[sample.name] = sample


    def _get_ariba_summary_data(self):
        ariba_summary_prefix = self.outprefix + '.ariba.summary'
        pyfastaq.utils.syscall('ariba summary --no_tree --preset cluster_small --known_variants ' + ariba_summary_prefix + ' ' + self.rootdir + '/*/ariba.report.tsv ')
        header = None
        self.ariba_raw_tsv_data = {}

        for line in open(ariba_summary_prefix + '.csv'):
            if header is None:
                header = line.rstrip().split(',')
            else:
                fields = line.rstrip().split(',')
                assert len(header) == len(fields)
                d = {x[0]: x[1] for x in zip(header[1:], fields[1:])}
                self.ariba_raw_tsv_data[fields[0]] = d


    def _load_ariba_results(self):
        self._get_ariba_summary_data()

        for sample in self.samples:
            sample_filename = os.path.join(self.rootdir, sample, 'ariba.report.tsv')
            self.samples[sample].tool_data['ariba'] = self.ariba_raw_tsv_data[sample_filename]
            self.samples[sample].load_ariba_report(sample_filename)


    def _load_kmerres_results(self):
        for sample in self.samples:
            filename = os.path.join(self.rootdir, sample, 'kmerResistance.out1')
            self.samples[sample].load_kmerres_file(filename)



    def _load_srst2_results(self):
        for sample in self.samples:
            fullgenes_file = os.path.join(self.rootdir, sample, 'srst2.out__fullgenes__results.txt')
            genes_file = os.path.join(self.rootdir, sample, 'srst2.out__genes__results.txt')
            self.samples[sample].load_srst2_genes_file(fullgenes_file, genes_file)


    def _load_all_data(self):
        self._load_ref_tsv()
        self._load_resistance_genes()
        self._load_ariba_results()
        self._load_kmerres_results()
        self._load_srst2_results()


    def _call_resistance(self):
        for sample in self.samples.values():
            sample.call_resistance(self.resistance_genes)

        self.antibio_calls = {}
        antibios = ['Amp', 'Cmp', 'GyrA', 'Str', 'Sul', 'Tet', 'Tri']
        self.venn = {}
        depths_of_amp_calls_only_by_kmerres = []
        kmerres_depths = []

        f = open(self.outprefix + '.calls.tsv', 'w')
        print('Antibiotic', 'Sample', 'Phenotype', '\t'.join(self.tools), 'ARIBA read depth', 'kmerResistance read depth', 'SRST2 read depth', sep='\t', file=f)

        for antibio in antibios:
            this_antibio_calls = {tool: {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0} for tool in self.tools}
            this_antibio_calls['R'] = 0
            this_antibio_calls['S'] = 0

            for sample_name, sample in self.samples.items():
                has_res = sample.ref_data[antibio]
                calls = {tool: sample.antibio_calls[tool][antibio] for tool in self.tools if tool != 'Ref'}
                calls['Ref'] = sample.ref_genotypes[antibio]
                tsv_data = ['1' if calls[tool] else '0' for tool in self.tools]

                if antibio == 'GyrA':
                    has_res ='S' if has_res == 'no' else 'R'
                else:
                    assert has_res in {'R', 'S', 'ND'}

                ariba_depth = 'NA' if sample.antibio_depths['ariba'][antibio] is None else round(sample.antibio_depths['ariba'][antibio], 1)
                kmerres_depth = 'NA' if sample.antibio_depths['kmerResistance'][antibio] is None else round(sample.antibio_depths['kmerResistance'][antibio], 1)
                srst2_depth = 'NA' if sample.antibio_depths['srst2'][antibio] is None else round(sample.antibio_depths['srst2'][antibio], 1)
                print(antibio, sample_name, has_res, '\t'.join(tsv_data), ariba_depth, kmerres_depth, srst2_depth, sep='\t', file=f)

                venn_key = tuple([calls[tool] for tool in self.tools])

                if venn_key not in self.venn:
                    self.venn[venn_key] = set()
                self.venn[venn_key].add(antibio + '.' + sample_name)


                if has_res == 'ND':
                    continue
                else:
                    this_antibio_calls[has_res] += 1


                for tool in calls:
                    if has_res == 'R':
                        if calls[tool]:
                            this_antibio_calls[tool]['TP'] += 1
                        else:
                            this_antibio_calls[tool]['FN'] += 1
                    else:
                        assert has_res == 'S'
                        if calls[tool]:
                            this_antibio_calls[tool]['FP'] += 1
                        else:
                            this_antibio_calls[tool]['TN'] += 1


                if antibio == 'Amp' and calls['kmerResistance'] and (not calls['ariba']) and (not calls['srst2']) and (not calls['Ref']):
                    depths_of_amp_calls_only_by_kmerres.append(sample.antibio_depths['kmerResistance']['Amp'])

                if sample.antibio_depths['kmerResistance'][antibio] is not None:
                    kmerres_depths.append(sample.antibio_depths['kmerResistance'][antibio])


            self.antibio_calls[antibio] = this_antibio_calls

        f.close()

        depths_of_amp_calls_only_by_kmerres.sort()
        print('Depth of Amp calls made only kmerResistance:', end=' ')
        print(*depths_of_amp_calls_only_by_kmerres, sep=', ')


    def _write_venn_file(self, filename):
        with open(filename, 'w') as f:
            print('\t'.join(self.tools), 'Count', 'Samples', sep='\t', file=f)

            for key, samples in self.venn.items():
                tool_columns = [{False: '0', True: '1'}[x] for x in key]
                print('\t'.join(tool_columns), len(samples), ';'.join(sorted(list(samples))), sep='\t', file=f)


    def _write_raw_results_to_file(self, filename):
        with open(filename, 'w') as f:
            for sample_name in sorted(self.samples):
                print(self.samples[sample_name], file=f)


    @staticmethod
    def _svg_circle(x, y, r, fill, stroke, stroke_width=1):
        return '<circle fill="' + fill + '" ' \
                + 'stroke="' + stroke + '" ' \
                + 'stroke-width="' + str(stroke_width) + '" ' \
                + 'cx="' + str(x) + '" ' \
                + 'cy="' + str(y) + '" ' \
                + 'r="' + str(r) + '" ' \
                + '/>'


    @staticmethod
    def _svg_line(x1, y1, x2, y2, color, thickness):
        return ''.join([
            '<line ',
            'stroke="' + color + '" ',
            'stroke-width="' + str(thickness) + '" ',
            'x1="' + str(x1) + '" ',
            'y1="' + str(y1) + '" ',
            'x2="' + str(x2) + '" ',
            'y2="' + str(y2) + '" ',
            '/>'])

    @staticmethod
    # coords  = list of tuples [(x1, y1), (x2, y2) ...]
    def _svg_polygon(coords, fill_colour, border_colour, border_width=1, opacity=-1):
        return_string = '<polygon points="' + ' '.join([str(x[0])+','+str(x[1]) for x in coords]) + '" ' \
                + 'fill="' + fill_colour + '" '

        if opacity != -1:
            return_string += 'fill-opacity="' + str(opacity) + '" '

        return_string += 'stroke="' + border_colour + '" ' \
                + 'stroke-width="' + str(border_width) + '" ' \
                + '/>'

        return return_string


    @staticmethod
    def _svg_text(x, y, text, fontsize, position='middle', writing_mode='lr', vertical=False, font_family='arial'):
        if vertical:
            vert = 'transform="rotate(-90,' + str(x) + ',' + str(y) + ')" '
        else:
            vert = ''

        return ''.join([
            '<text ',
            vert,
            'text-anchor="' + position + '" ',
            'font-size="' + str(fontsize) + '" ',
            'font-family="' + font_family + '" ',
            'writing-mode="' + writing_mode + '" ',
            'x="' + str(x) + '" ',
            'y="' + str(y) + '">',
            text,
            '</text>'])


    def _make_upset_plot(self):
        f = open(self.outprefix + '.upset.svg', 'w')
        width = 650
        height = 600
        plot_left_start = 150
        plot_right_end = 590
        bar_plot_y_top = 30
        bar_plot_origin_x = plot_left_start
        bar_plot_origin_y = 450
        bar_x_gap = 10

        print(r'''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg width="''' + str(width) + '" height="' + str(height) + '">', file=f)

        venn_sizes = {x: len(self.venn[x]) for x in self.venn}

        number_of_bars = len(venn_sizes)
        bar_x_step = (plot_right_end - plot_left_start) / number_of_bars
        bar_width = bar_x_step - bar_x_gap
        bar_y_max = 550
        number_of_dot_rows = len(self.tools)
        dot_radius = 0.5 * bar_width
        if self.lenient:
            dot_radius *= 0.75
        top_dot_y = bar_plot_origin_y + (2 * dot_radius)
        bottom_dot_y = height - (2 * dot_radius)
        dot_y_step = (bottom_dot_y - top_dot_y) / (number_of_dot_rows - 1)
        x_position = bar_plot_origin_x

        # --------------  Bars and circles ---------------#
        for key, set_size in sorted(venn_sizes.items(), key=operator.itemgetter(1), reverse=True):
            # ---------- Bars --------#
            bar_x_start = x_position + bar_x_gap
            bar_x_end = bar_x_start + bar_width
            bar_x_centre = bar_x_start + 0.5 * (bar_x_end - bar_x_start)
            bar_height = (bar_plot_origin_y - bar_plot_y_top) * set_size / bar_y_max
            bar_y_top = bar_plot_origin_y - bar_height
            bar_coords = [(bar_x_start, bar_plot_origin_y), (bar_x_start, bar_y_top), (bar_x_end, bar_y_top), (bar_x_end, bar_plot_origin_y)]
            print(self._svg_polygon(bar_coords, "darkgray", "black"), file=f)
            print(self._svg_text(bar_x_centre, bar_y_top - 5, str(set_size), 12, position="middle"), file=f)

            # ---------- Circles --------#
            y_position = top_dot_y
            for tool, present in zip(self.tools, key):
                if present:
                    fill_colour = tool_colours[tool]
                    border_colour = "black"
                else:
                    fill_colour = "lightgray"
                    border_colour = "lightgray"
                print(self._svg_circle(bar_x_centre, y_position, dot_radius, fill_colour, border_colour, stroke_width=3), file=f)
                y_position += dot_y_step

            x_position += bar_x_step

        # --------------  y axis ---------------#
        print(self._svg_text(bar_plot_origin_x - 40, 0.5 * (bar_plot_y_top + bar_plot_origin_y), 'Size of intersection', 14, position='middle', vertical=True), file=f)
        print(self._svg_line(bar_plot_origin_x, bar_plot_origin_y, bar_plot_origin_x, bar_plot_y_top, "black", 2), file=f)
        for i in range(0, 501, 100):
            y = bar_plot_origin_y - (bar_plot_origin_y - bar_plot_y_top) * i / bar_y_max
            print(self._svg_line(bar_plot_origin_x - 3, y, bar_plot_origin_x, y, "black", 1), file=f)
            print(self._svg_text(bar_plot_origin_x - 7, y+3, str(i), 10, position='end'), file=f)

        # --------------  tool labels ---------------#
        y = top_dot_y
        tool_translate = {'ariba':'ARIBA', 'srst2':'SRST2', 'kmerResistance': 'KmerResistance', 'Ref': 'Holt 2012'}
        for tool in self.tools:
            print(self._svg_text(bar_plot_origin_x - 5, y+6, tool_translate[tool], 16, position='end'), file=f)
            y += dot_y_step


        print('</svg>', file=f)

        f.close()


    def run(self):
        self._load_all_data()
        self._write_amr_genes_suppl_tsv()
        self._call_resistance()
        self._write_raw_results_to_file(self.outprefix + '.debug.raw_results.tsv')
        self._write_venn_file(self.outprefix + '.venn.tsv')
        self._make_upset_plot()


parser = argparse.ArgumentParser(
    description = 'Makes plots from shigella data',
)

parser.add_argument('--lenient', action='store_true', help='Also count SRST2 with ? as present, and ariba partial/fragmented/interrupted as present')
parser.add_argument('rootdir', help='Shigella/ directory')
parser.add_argument('srst2_ref_fa', help='SRST2 reference fasta')
parser.add_argument('ariba_db_dir', help='ARIBA db directory')
parser.add_argument('outprefix', help='Prefix of output files')
options = parser.parse_args()


ref_data_tsv = os.path.join(options.rootdir, 'shigella_supplementary_holt2012.tsv')
ariba_clusters_tsv = os.path.join(options.ariba_db_dir, '02.cdhit.clusters.tsv')
samples = Samples(ref_data_tsv, options.srst2_ref_fa, ariba_clusters_tsv, options.rootdir, options.outprefix, options.lenient)
samples.run()
