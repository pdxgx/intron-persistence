#!/usr/bin/env python3

"""
intronomer.py
Python 3 code to calculate intron persistence for a long-read RNA-seq sample

Step 1 assigns aligned pacbio reads to annotated transcripts
Step 2 calculates intron persistence values from read-transcript sets

SAMPLE RUN:
intron-persistence/intronomer/intronomer.py
    -g <annotation file>
    -a <aligned-reads file>
    -p <project flag>
    -o <output directory>
    -t <read-transcript preprocess file>

<aligned-reads file> (.bam or .sam format)
    - output of long RNA-seq reads mapped to a reference with a splice-aware
      aligner
    - If using PacBio reads, full length non-concatamer poly-a selected long
      reads are recommended; please do not use reads that have been clustered
      and collapsed into a consensus isoform.
<annotation file> (.gtf format)
    - annotation based on the same reference genome to which the long RNA-seq
      reads were mapped. Any potentially retained persistent introns will be
      defined with respect to this annotation.
<project flag> (optional)
    - a flag or project name to add to output filenames
<output directory> (optional)
    - desired output directory; if not given, the alignment file directory
      will be used to write new files
<read-transcript preprocess file> (optional)
    - if the RI_txs_to_read_ids file was previously generated (first step in
      the calculation) you can add it here to avoid re-generating this.
"""
import argparse
from collections import defaultdict
from datetime import datetime
from math import ceil
import numpy as np
import os
import pandas as pd
import pysam
from scipy.spatial.distance import cdist


_STRAND_MAP = {True: '-', False: '+'}


def create_intron_matrix(transcript_df):
    transcript_df.sort_values(
        by="read_introns", key=lambda x: x.str.len(),
        ascending=False, inplace=True
    )
    tx_introns = transcript_df['tx_introns'].unique()[0]
    tx_introns = [intron for intron in tx_introns.split(';') if intron]
    if transcript_df['strand'].unique()[0] == '-':
        tx_introns = tx_introns[::-1]
    plot_list = []
    reads = []
    for index, row in transcript_df.iterrows():
        plot_introns = []
        reads.append(row['read_id'])
        read_ints = row['read_introns'].split(';')
        for intron in tx_introns:
            if intron in read_ints:
                plotval = 1
            else:
                i_left, i_right = map(int, intron.split('-'))
                intron_in_read = (
                        i_left > row['read_left']
                        and i_right < row['read_right']
                )
                if intron_in_read:
                    plotval = 0
                else:
                    plotval = np.NaN
            plot_introns.append(plotval)
        plot_list.append(plot_introns)

    intron_df = pd.DataFrame(
        np.array(plot_list), columns=tx_introns, index=reads
    )
    return intron_df


def hamming_similarity(df_row, df):
    hs = []
    df_row.rename('orig', inplace=True)
    for index, row in df.iterrows():
        r_df = df_row.to_frame().join(row)
        r_df.dropna(axis=0, inplace=True)
        hs.append(
            cdist(np.array([r_df['orig']]), np.array([r_df[index]]), 'hamming')
        )
    return 1 - np.mean(hs)


def assign_RI_metrics(read_tx_df, output_dir, now, batch_num, flag=''):
    read_tx_df = read_tx_df.loc[read_tx_df['RI-read_ratio'] > 0].copy()
    all_intron_columns = [
        'intron', 'transcript', 'persistence', 'position', 'chrom'
    ]
    all_intron_rows = []
    for tx in read_tx_df['transcript'].unique():
        sub_df = read_tx_df.loc[read_tx_df['transcript'] == tx].copy()
        chrom = sub_df['chrom'].unique()[0]
        intron_df = create_intron_matrix(sub_df)
        num_introns = len(intron_df.columns) - 1
        n_x = len(intron_df)
        for int_count, intron in enumerate(intron_df.columns.values):
            int_df = intron_df.dropna(subset=[intron]).copy()
            if len(int_df) == 0:
                persistence = 0
            elif int_df[intron].sum() == int_df[intron].count():
                persistence = 0
            else:
                # Percent of transcript reads with intron information
                info_perc = len(int_df) / n_x
                # Only need to calculate splicing progression and similarity
                # if the intron is retained in the read
                int_df = int_df.loc[int_df[intron] == 0]
                # Splicing progression ignoring the target intron
                int_df['p_xy'] = (
                    int_df.drop([intron], axis=1).sum(axis=1)
                    / int_df.drop([intron], axis=1).count(axis=1)
                )
                # Setup for multiplying in read-specific hamming distance
                int_df['h_sim'] = int_df.apply(
                    lambda x: hamming_similarity(x, int_df), axis=1
                )
                int_df['IR_persistence'] = int_df.apply(
                    lambda x: x['p_xy'] * (1 - x[intron]) * x['h_sim'],
                    axis=1
                )
                persistence = info_perc * int_df['IR_persistence'].sum() / n_x
            all_intron_rows.append((
                f'{chrom}:{intron}', tx, persistence, int_count / num_introns,
                chrom
            ))

    tx_df = pd.DataFrame(columns=all_intron_columns, data=all_intron_rows)
    print_df = tx_df[['intron', 'transcript', 'persistence', 'position']]
    out_name = f'intron-persistence_{now}.csv'
    if flag:
        out_name = f'{flag}_{out_name}'
    print_header=True
    if batch_num is not None:
        if batch_num > 0:
            print_header = False
        out_name = f'batch{batch_num}_' + out_name
    outfile = os.path.join(output_dir, out_name)
    print_df.to_csv(outfile, index=False, sep=',', header=print_header)
    return


def extract_introns(gtf_file):
    """Extracts splice site annotations from .gtf file
    This function is a modified version of one that is part of HISAT.
    Copyright 2014, Daehwan Kim <infphilo@gmail.com>
    HISAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    HISAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
    """
    genes = defaultdict(list)
    trans = {}
    intron_to_txs = {}
    tx_ranges = {}
    tx_to_introns = {}
    tx_to_gene = {}
    with open(gtf_file) as gtf:
        # Parse valid exon lines from the annotation file into a dict by
        # transcript_id
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()

            try:
                items = line.split('\t')
                chrom, _, feature, left, right, _, strand, _, values = items
            except ValueError:
                continue
            left, right = int(left) - 1, int(right) + 1

            if feature != 'exon' or left >= right:
                continue

            vals_dict = {}
            for attr in values.split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                vals_dict[attr] = val.strip('"')

            if 'gene_id' not in vals_dict or 'transcript_id' not in vals_dict:
                continue

            transcript_id = vals_dict['transcript_id']
            if transcript_id not in trans:
                trans[transcript_id] = [chrom, strand, [[left, right]]]
                genes[vals_dict['gene_id']].append(transcript_id)
            else:
                trans[transcript_id][2].append([left, right])
            tx_to_gene[transcript_id] = vals_dict['gene_id']

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions
    for tx_name, (chrom, strand, exons) in trans.items():
        tx_to_introns[tx_name] = set()
        if chrom not in intron_to_txs:
            intron_to_txs[chrom] = {}
            intron_to_txs[chrom]['+'] = defaultdict(set)
            intron_to_txs[chrom]['-'] = defaultdict(set)
            intron_to_txs[chrom]['.'] = defaultdict(set)
        for i in range(1, len(exons)):
            tx_to_introns[tx_name].add((exons[i - 1][1], exons[i][0]))
            left = exons[i - 1][1]
            right = exons[i][0]
            intron_to_txs[chrom][strand][(left, right)].add(tx_name)
        tx_ranges[tx_name] = [exons[0][0], exons[-1][1]]
        tx_to_introns[tx_name] = sorted(tx_to_introns[tx_name])

    return tx_ranges, intron_to_txs, tx_to_introns, tx_to_gene


def check_match_type(read_introns, read_left, read_right, tx_introns):
    splice_type = 'all_introns-'
    read_length = 'full_length'
    target_introns = 0
    missing_introns = 0
    for i, intron in enumerate(tx_introns, 1):
        if intron in read_introns:
            target_introns += 1
            continue
        if intron[0] > read_left and intron[1] < read_right:
            splice_type = 'skipped_splicing-'
            target_introns += 1
            missing_introns += 1
        else:
            read_length = 'partial'
    return splice_type + read_length, target_introns, missing_introns


def longreads_to_isoforms(bam_file, intron_info, output_dir, now, flag):
    tx_ranges, intron_to_txs, tx_to_introns, tx_to_gene = intron_info
    if bam_file.endswith('bam'):
        # print('bamfile')
        samfile = pysam.AlignmentFile(bam_file, 'rb')
    else:
        # print('samfile')
        samfile = pysam.AlignmentFile(bam_file, 'r')
    read_maps_encountered = set()
    read_info_rows = []
    read_info_cols = [
        'read_id', 'orig_read_id', 'read_length', 'read_left', 'read_right',
        'chrom', 'strand', 'read_introns',
    ]
    tx_info_rows = []
    tx_info_cols = [
        'transcript', 'tx_length', 'tx_left', 'tx_right', 'tx_introns'
    ]
    tx_to_read_id_rows = []
    tx_to_read_id_cols = [
        'transcript', 'read_id', 'match_type', 'region_tx_intron_count',
        'unspliced_intron_count'
    ]
    fl_match = 'all_introns-full_length'
    for i, read in enumerate(samfile):
        collect_read = False
        read_str = read.to_string()
        if read_str in read_maps_encountered:
            continue
        assigned_read_number = i
        read_maps_encountered.add(read_str)
        possible_txs = set()
        cigar_tups = read.cigartuples
        if not cigar_tups:
            continue
        left_pos = read.reference_start
        right_pos = read.reference_end
        strand = _STRAND_MAP[read.is_reverse]
        chrom = read.reference_name
        curr_pos = left_pos
        introns = []
        for (c_type, length) in cigar_tups:
            if c_type == 3:
                i_left = curr_pos + 1
                curr_pos += length
                i_right = curr_pos
                introns.append((i_left, i_right))
            elif c_type in (0, 2):
                curr_pos += length
        for intron in introns:
            try:
                possible_txs.update(intron_to_txs[chrom][strand][intron])
            except KeyError:
                continue
        intron_set = set(introns)
        r_intron_str = ''
        for (left, right) in sorted(list(intron_set)):
            r_intron_str += f'{left}-{right};'
        shortlist_txs = set()
        exact_txs = set()
        for tx in possible_txs:
            target_intron_set = set(tx_to_introns[tx])
            if intron_set == target_intron_set:
                exact_txs.add(tx)
                collect_read = True

        exact_tx_count = len(exact_txs)
        if exact_txs:
            if exact_tx_count == 1:
                best_tx = list(exact_txs)[0]
            else:
                min_len_diff = 1000000000000
                best_tx = ''
                for tx in exact_txs:
                    tx_left = tx_ranges[tx][0]
                    tx_right = tx_ranges[tx][1]
                    len_diff = abs(
                        (tx_right - tx_left) - (right_pos - left_pos)
                    )
                    if len_diff < min_len_diff:
                        best_tx = tx
                        min_len_diff = len_diff
            tx_to_read_id_rows.append((
                best_tx, assigned_read_number, fl_match,
                len(set(tx_to_introns[best_tx])), 0
            ))
        else:
            for tx in possible_txs:
                target_intron_set = set(tx_to_introns[tx])
                if intron_set.issubset(target_intron_set):
                    shortlist_txs.add(tx)
                    collect_read = True

            if shortlist_txs:
                if len(shortlist_txs) == 1:
                    best_tx = list(shortlist_txs)[0]
                else:
                    min_len_diff = 1000000000000
                    best_tx = ''
                    for tx in shortlist_txs:
                        tx_left = tx_ranges[tx][0]
                        tx_right = tx_ranges[tx][1]
                        len_diff = abs(
                            (tx_right - tx_left) - (right_pos - left_pos)
                        )
                        if len_diff < min_len_diff:
                            best_tx = tx
                            min_len_diff = len_diff
                match_info = check_match_type(
                    sorted(list(intron_set)), left_pos, right_pos,
                    sorted(list(tx_to_introns[best_tx]))
                )
                match_type, tx_count, missing_introns = match_info
                tx_to_read_id_rows.append((
                    best_tx, assigned_read_number, match_type, tx_count,
                    missing_introns
                ))
        if collect_read:
            read_info_rows.append((
                assigned_read_number, read.query_name, read.reference_length,
                left_pos, right_pos, chrom, strand, r_intron_str
            ))

    read_info_df = pd.DataFrame(
        columns=read_info_cols, data=read_info_rows
    )
    # outfile = os.path.join(
    #     output_dir, f'read_info_nosequence_{flag}_{now}.tsv'
    # )
    # with open(outfile, 'w') as output:
    #     read_info_df.to_csv(output, sep='\t', index=False)

    tx_info_df = pd.DataFrame(
        columns=tx_info_cols, data=tx_info_rows
    )
    # outfile = os.path.join(
    #     output_dir, f'tx_info_{flag}_{now}.tsv'
    # )
    # with open(outfile, 'w') as output:
    #     tx_info_df.to_csv(output, sep='\t', index=False)

    for tx in tx_info_df['transcript'].unique():
        target_intron_set = set(tx_to_introns[tx])
        t_intron_str = ''
        for (left, right) in sorted(list(target_intron_set)):
            t_intron_str += f'{left}-{right};'
        tx_info_rows.append((
            tx, tx_ranges[tx][1] - tx_ranges[tx][0], tx_ranges[tx][0],
            tx_ranges[tx][1], t_intron_str
        ))
    tx_to_read_id_df = pd.DataFrame(
        columns=tx_to_read_id_cols, data=tx_to_read_id_rows
    )
    # outfile = os.path.join(
    #     output_dir, f'tx_to_read_id_{flag}_{now}.tsv'
    # )
    # tx_to_read_id_df['gene_id'] = tx_to_read_id_df['transcript'].apply(
    #     lambda x: tx_to_gene.get(x, '')
    # )
    # with open(outfile, 'w') as output:
    #     tx_to_read_id_df.to_csv(output, sep='\t', index=False)

    # tx_groups = tx_to_read_id_df.groupby('transcript').count()
    # tx_groups['reads_per_transcript'] = tx_groups['read_id']
    # tx_groups.reset_index(inplace=True)
    # tx_to_gene = tx_to_read_id_df.set_index('transcript')['gene_id'].to_dict()
    # tx_groups['gene_id'] = tx_groups['transcript'].apply(
    #     lambda x: tx_to_gene[x]
    # )
    # tx_groups = tx_groups[['gene_id', 'transcript', 'reads_per_transcript']]
    # gene_groups = tx_to_read_id_df.groupby('gene_id').count()
    # gene_groups['reads_per_gene'] = gene_groups['read_id']
    # gene_groups = gene_groups[['reads_per_gene']]
    # gene_groups.reset_index(inplace=True)
    # mergedf = pd.merge(gene_groups, tx_groups, how='outer', on='gene_id')
    # outfile = os.path.join(output_dir, 'reads_per_gene_and_transcript.tsv')
    # with open(outfile, 'w') as output:
    #     mergedf.to_csv(output, sep='\t', index=False)

    tx_to_read_id_df = pd.merge(
        tx_to_read_id_df, read_info_df, on=['read_id'], how='left'
    )
    tx_to_read_id_df = pd.merge(
        tx_to_read_id_df, tx_info_df, on=['transcript'], how='left'
    )
    ratios_dict = {}
    RI_tx_set = set()
    for tx in tx_to_read_id_df['transcript'].unique():
        mini_df = tx_to_read_id_df.loc[tx_to_read_id_df['transcript'] == tx]
        skips = len(
            mini_df.loc[mini_df['match_type'].isin(
                ['skipped_splicing-full_length', 'skipped_splicing-partial'])]
        )
        if skips > 0 and len(mini_df) > 4:
            RI_tx_set.add(tx)
            ratios_dict[tx] = skips / len(mini_df)
    tx_to_read_id_df = tx_to_read_id_df.loc[
        tx_to_read_id_df['transcript'].isin(RI_tx_set)
    ]
    tx_to_read_id_df['RI-read_ratio'] = tx_to_read_id_df['transcript'].apply(
        lambda x: f'{ratios_dict[x]:.2f}'
    )
    tx_to_read_id_df.sort_values(
        by=['match_type', 'transcript', 'RI-read_ratio'],
        ascending=[False, True, False], inplace=True
    )
    outfile = os.path.join(
        output_dir, f'RI_txs_to_read_ids_{flag}_{now}.tsv'
    )
    with open(outfile, 'w') as output:
        tx_to_read_id_df.to_csv(output, sep='\t', index=False)
    num_txs = tx_to_read_id_df['transcript'].nunique()
    print(
        f'\n{num_txs} transcripts with at least 5 reads and some persistent '
        f'RIs'
    )
    return tx_to_read_id_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculates intron persistence in long-read RNA-seq data.'
    )
    parser.add_argument(
        '--aligned-long-reads', '-a', required=True,
        help='Sorted .bam or .sam file for aligned long reads. If using '
             'PacBio reads, full length non-concatamer poly-a selected long '
             'reads are recommended; please do not use reads that have been '
             'clustered and collapsed into a consensus isoform.'
    )
    parser.add_argument(
        '--annotation-gtf', '-g', required=True,
        help='.gtf file with GENCODE annotation.'
    )
    parser.add_argument(
        '--batch-number', type=int,
        help='For running in batches: the current batch number (e.g. a slurm '
             'array task id) to process correct chunk of the dataframe.'
    )
    parser.add_argument(
        '--total-batches', type=int,
        help='If this is a batch run: this is the total number of batches'
    )
    parser.add_argument(
        '--project-flag', '-p',
        help='If desired, add a project name or other flag that you would '
             'like to have included in the output filenames.'
    )
    parser.add_argument(
        '--transcript-read-file', '-t',
        help='If the "RI_txs_to_read_ids" file has already been generated, '
             'skip re-generation by passing it here.'
    )
    parser.add_argument(
        '--output-directory', '-o',
        help='Add the desired directory in which to write output files. If '
             'none is provided, files will be written in the same directory '
             'as the input alignment file.'
    )

    args = parser.parse_args()
    bam_path = args.aligned_long_reads
    gtf_path = args.annotation_gtf
    batch_num = args.batch_number
    tot_batches = args.total_batches
    project_flag = args.project_flag
    read_tx_file = args.transcript_read_file
    output_dir = args.output_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    if output_dir is None:
        output_dir = os.path.dirname(bam_path)
    os.makedirs(output_dir, exist_ok=True)
    if read_tx_file is None:
        intron_info = extract_introns(gtf_path)
        read_tx_df = longreads_to_isoforms(
            bam_path, intron_info, output_dir, now, project_flag
        )
    else:
        read_tx_df = pd.read_csv(read_tx_file, sep='\t')
    if batch_num is not None:
        all_txs = read_tx_df['transcript'].unique().tolist()
        all_txs.sort()
        chunksize = ceil(len(all_txs) / tot_batches)
        tx_lists = [
            all_txs[i:i + chunksize]
            for i in range(0, len(all_txs), chunksize)
        ]
        target_txs = tx_lists[batch_num]
        print(
            f'Batch {batch_num}: {len(all_txs)} total transcripts, '
            f'{len(target_txs)} in current batch'
        )
        read_tx_df = read_tx_df.loc[read_tx_df['transcript'].isin(target_txs)]
        print(
            f"total length of df: {len(read_tx_df)}, "
            f"{read_tx_df['transcript'].nunique()} txs"
        )
    assign_RI_metrics(read_tx_df, output_dir, now, batch_num, project_flag)
