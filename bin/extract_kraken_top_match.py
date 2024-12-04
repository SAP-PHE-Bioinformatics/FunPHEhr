#!/usr/bin/env python
import os
import sys
import argparse
import gzip
from time import gmtime, strftime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Tree:
    'Tree node.'
    def __init__(self, taxid, level_num, level_id, children=None, parent=None):
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)

def process_kraken_output(kraken_line):
    l_vals = kraken_line.split('\t')
    if len(l_vals) < 5:
        return [-1, '']
    tax_id = int(l_vals[2].split('taxid ')[-1][:-1] if "taxid" in l_vals[2] else l_vals[2])
    read_id = l_vals[1]
    return [tax_id, read_id]

def process_kraken_report(report_line):
    l_vals = report_line.strip().split('\t')
    if len(l_vals) < 5:
        return []
    try:
        int(l_vals[1])
    except ValueError:
        return []
    try:
        taxid = int(l_vals[-3])
        level_type = l_vals[-2]
        map_kuniq = {'species': 'S', 'genus': 'G', 'family': 'F',
                     'order': 'O', 'class': 'C', 'phylum': 'P',
                     'superkingdom': 'D', 'kingdom': 'K'}
        if level_type not in map_kuniq:
            level_type = '-'
        else:
            level_type = map_kuniq[level_type]
    except ValueError:
        taxid = int(l_vals[-2])
        level_type = l_vals[-3]
    spaces = len(l_vals[-1]) - len(l_vals[-1].lstrip())
    level_num = int(spaces / 2)
    return [taxid, level_num, level_type]

def parse_report(report_file):
    taxids = {}
    tree_nodes = {}
    prev_node = None
    with open(report_file, 'r') as r_file:
        for line in r_file:
            report_vals = process_kraken_report(line)
            if len(report_vals) == 0:
                continue
            taxid, level_num, level_id = report_vals
            if taxid == 0 or taxid == 1:
                continue
            node = Tree(taxid, level_num, level_id)
            if prev_node is None:
                prev_node = node
                tree_nodes[taxid] = node
                continue
            while prev_node and prev_node.level_num >= node.level_num:
                prev_node = prev_node.parent
            if prev_node:
                prev_node.add_child(node)
                node.parent = prev_node
            tree_nodes[taxid] = node
            prev_node = node
            taxids[taxid] = node
    return taxids, tree_nodes

def find_top_family(taxids, report_file):
    family_taxid = None
    family_count = 0
    with open(report_file, 'r') as r_file:
        for line in r_file:
            report_vals = process_kraken_report(line)
            if len(report_vals) == 0:
                continue
            taxid, level_num, level_id = report_vals
            if level_id == 'F':
                count = int(line.strip().split('\t')[1])
                if count > family_count:
                    family_count = count
                    family_taxid = taxid
    return family_taxid

def get_descendants(node):
    descendants = set()
    nodes_to_visit = [node]
    while nodes_to_visit:
        current_node = nodes_to_visit.pop()
        descendants.add(current_node.taxid)
        nodes_to_visit.extend(current_node.children)
    return descendants

def get_ancestors_up_to_family(node):
    ancestors = set()
    current_node = node
    while current_node.parent and current_node.parent.level_id != 'F':
        ancestors.add(current_node.parent.taxid)
        current_node = current_node.parent
    if current_node.parent and current_node.parent.level_id == 'F':
        ancestors.add(current_node.parent.taxid)
    return ancestors

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', dest='kraken_file', required=True,
        help='Kraken output file to parse')
    parser.add_argument('-s','-s1', '-1', '-U', dest='seq_file1', required=True,
        help='FASTA/FASTQ File containing the raw sequence letters.')
    parser.add_argument('-s2', '-2', dest='seq_file2', default= "",
        help='2nd FASTA/FASTQ File containing the raw sequence letters (paired).')
    parser.add_argument('-o', "--output", dest='output_file', required=True,
        help='Output FASTA/Q file containing the combined reads and sample IDs')
    parser.add_argument('--append', dest='append', action='store_true',
        help='Append the sequences to the end of the output FASTA file specified.')
    parser.add_argument('--noappend', dest='append', action='store_false',
        help='Create a new FASTA file containing sample sequences and IDs (rewrite if existing) [default].')
    parser.add_argument('--max', dest='max_reads', required=False, 
        default=100000000, type=int,
        help='Maximum number of reads to save [default: 100,000,000]')
    parser.add_argument('-r', '--report', dest='report_file', required=True,
        help='Kraken report file.')
    parser.add_argument('--fastq-output', dest='fastq_out', required=False,
        action='store_true', default=False,
        help='Print output FASTQ reads [requires input FASTQ, default: output is FASTA]')
    parser.set_defaults(append=False)

    args = parser.parse_args()
    
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')
    
    taxids, tree_nodes = parse_report(args.report_file)
    if not taxids:
        sys.stderr.write("No taxids found in report file\n")
        sys.exit(1)

    top_family_taxid = find_top_family(taxids, args.report_file)
    if top_family_taxid is None:
        sys.stderr.write("No family found in report file\n")
        sys.exit(1)

    save_taxids = set()
    descendants = get_descendants(tree_nodes[top_family_taxid])
    ancestors = get_ancestors_up_to_family(tree_nodes[top_family_taxid])
    save_taxids.update(descendants)
    save_taxids.update(ancestors)
    save_taxids.add(top_family_taxid)

    sys.stdout.write(f"\tTaxonomy IDs to parse: {', '.join(map(str, save_taxids))}\n")
    sys.stdout.write(">> STEP 1: PARSING KRAKEN FILE FOR READIDS %s\n" % args.kraken_file)
    count_kraken = 0
    read_line = -1
    save_readids = set()
    k_file = open(args.kraken_file, 'r')
    sys.stdout.write('\t0 reads processed')
    sys.stdout.flush()
    for line in k_file:
        count_kraken += 1
        if (count_kraken % 10000 == 0):
            sys.stdout.write('\r\t%0.2f million reads processed' % float(count_kraken/1000000.))
            sys.stdout.flush()
        tax_id, read_id = process_kraken_output(line)
        if tax_id == -1:
            continue
        if tax_id in save_taxids:
            save_readids.add(read_id)
    k_file.close()
    sys.stdout.write('\r\t%0.2f million reads processed\n' % float(count_kraken/1000000.))
    sys.stdout.write('\t%i read IDs saved\n' % len(save_readids))

    seq_file1 = args.seq_file1
    seq_file2 = args.seq_file2
    if(seq_file1[-3:] == '.gz'):
        s_file1 = gzip.open(seq_file1, 'rt')
    else:
        s_file1 = open(seq_file1, 'rt')
    first = s_file1.readline()
    if len(first) == 0:
        sys.stderr.write("ERROR: sequence file's first line is blank\n")
        sys.exit(1)
    if first[0] == ">":
        filetype = "fasta"
    elif first[0] == "@":
        filetype = "fastq"
    else:
        sys.stderr.write("ERROR: sequence file must be FASTA or FASTQ\n")
        sys.exit(1)
    s_file1.close()
    if filetype != 'fastq' and args.fastq_out:
        sys.stderr.write('ERROR: for FASTQ output, input file must be FASTQ\n')
        sys.exit(1)
    if(seq_file1[-3:] == '.gz'):
        s_file1 = gzip.open(seq_file1, 'rt')
        if len(seq_file2) > 0:
            s_file2 = gzip.open(seq_file2, 'rt')
    else:
        s_file1 = open(seq_file1, 'r')
        if len(seq_file2) > 0:
            s_file2 = open(seq_file2, 'r')

    sys.stdout.write(">> STEP 2: READING SEQUENCE FILES AND WRITING READS\n")
    sys.stdout.write('\t0 read IDs found (0 mill reads processed)')
    sys.stdout.flush()

    if args.append:
        o_file = open(args.output_file, 'a')
    else:
        o_file = open(args.output_file, 'w')

    count_seqs = 0
    count_output = 0
    for record in SeqIO.parse(s_file1, filetype):
        count_seqs += 1
        if (count_seqs % 1000 == 0):
            sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
            sys.stdout.flush()
        test_id = str(record.id)
        test_id2 = test_id
        if ("/1" in test_id) or ("/2" in test_id):
            test_id2 = test_id[:-2]
        if test_id in save_readids or test_id2 in save_readids:
            count_output += 1
            sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
            sys.stdout.flush()
            SeqIO.write(record, o_file, "fastq" if args.fastq_out else "fasta")
    s_file1.close()
    o_file.close()
    sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)\n' % (count_output, float(count_seqs/1000000.)))
    sys.stdout.flush()
    
    if len(seq_file2) > 0:
        if(seq_file2[-3:] == '.gz'):
            s_file2 = gzip.open(seq_file2, 'rt')
        else:
            s_file2 = open(seq_file2, 'r')
        if args.append:
            o_file = open(args.output_file, 'a')
        else:
            o_file = open(args.output_file, 'w')
        count_output = 0
        count_seqs = 0
        sys.stdout.write('\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
        sys.stdout.flush()
        for record in SeqIO.parse(s_file2, filetype):
            count_seqs += 1
            if (count_seqs % 1000 == 0):
                sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
                sys.stdout.flush()
            test_id = str(record.id)
            test_id2 = test_id
            if ("/1" in test_id) or ("/2" in test_id):
                test_id2 = test_id[:-2]
            if test_id in save_readids or test_id2 in save_readids:
                count_output += 1
                sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)' % (count_output, float(count_seqs/1000000.)))
                sys.stdout.flush()
                SeqIO.write(record, o_file, "fastq" if args.fastq_out else "fasta")
        s_file2.close()
        o_file.close()
        sys.stdout.write('\r\t%i read IDs found (%0.2f mill reads processed)\n' % (count_output, float(count_seqs/1000000.)))
    
    sys.stdout.write('\tGenerated file: %s\n' % args.output_file)
    
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')
    sys.exit(0)

if __name__ == "__main__":
    main()