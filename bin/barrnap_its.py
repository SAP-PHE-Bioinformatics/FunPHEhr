#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import subprocess
import sys
import pandas as pd

def get_params(argv):
    parser = argparse.ArgumentParser(description='Extract ITS1 from fungal genome rapidly')
    parser.add_argument('-i', '--i', help="input genome file", required=True)
    parser.add_argument('-o', '--o', help="output directory", required=True)
    parser.add_argument('-which', '--which', help="Which ITS sequence to extract (ITS1|ITS2) default=ITS1", default='ITS1')
    parser.add_argument('-cpu', '--cpu', help="number of threads/cores to use", required=False, default='48')
    parser.add_argument('-name', '--name', help="name", required=False, default='genome')
    a = parser.parse_args()
    return a

def terminal(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    p.wait()

if __name__ == '__main__':
    a = get_params(sys.argv[1:])
    if a.o[-1] != '/':
        a.o = a.o + '/'

    if a.name == 'genome':
        name = a.i.split('/')[-1].split('.')[0]
    else:
        name = a.name

    if a.which not in ['ITS1', 'ITS2', 'all']:
        print("ERROR, argument --which not recognized. Please choose between 'ITS1' and 'ITS2'")

    ## Running barrnap
    terminal('barrnap --kingdom euk --threads ' + a.cpu + ' ' + a.i + ' > ' + a.o + 'barrnap.tmp')

    ## Parsing and getting rDNA gene cluster coordinates
    gff = pd.read_csv(a.o + 'barrnap.tmp', header=None, comment='#', sep='\t')
    gff = gff[gff[8].str.contains('18S') | gff[8].str.contains('28S')].sort_values(by=[0, 6, 3, 4]).reset_index(drop=True)
    regions = []
    print(gff)
    for i in gff.index:
        try:
            if 'partial' in gff.loc[i, 8]:
                pass
            elif gff.loc[i, 6] == '+':
                if '18S_rRNA' in gff.loc[i, 8] and '28S_rRNA' in gff.loc[i + 1, 8] and gff.loc[i, 0] == gff.loc[i + 1, 0] and gff.loc[i, 6] == gff.loc[i + 1, 6]:
                    regions.append((gff.loc[i, 0], gff.loc[i, 6], gff.loc[i, 3] - 1, gff.loc[i + 1, 4] - 1))
            elif gff.loc[i, 6] == '-':
                if '28S_rRNA' in gff.loc[i, 8] and '18S_rRNA' in gff.loc[i + 1, 8] and gff.loc[i, 0] == gff.loc[i + 1, 0] and gff.loc[i, 6] == gff.loc[i + 1, 6]:
                    regions.append((gff.loc[i, 0], gff.loc[i, 6], gff.loc[i, 3] - 1, gff.loc[i + 1, 4] - 1))
            else:
                print('Error, impossible to read ' + gff.loc[i, 6])
        except:
            pass

    print(regions)

    ## Opening genome fasta, extracting regions
    with open(a.i, 'r') as handle:
        fasta = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))

    extractedRegions = []
    for r in regions:
        print(r)
        SEQ = fasta[r[0]][r[2]:r[3]]
        SEQ.id = r[0] + '_' + r[1] + '_' + str(r[2]) + '_' + str(r[3])
        SEQ.name = SEQ.id
        SEQ.description = SEQ.id
        extractedRegions.append(SEQ)
        print(SEQ)

    with open(a.o + 'regions.fasta', "w") as output_handle:
        SeqIO.write(extractedRegions, output_handle, 'fasta')

    if len(extractedRegions) > 0:
        ## Run ITSx on regions.fasta
        terminal('ITSx -t F -i %s -o %s --cpu %s --save_regions %s --not_found F --graphical F --summary F  --fasta F' % (a.o + 'regions.fasta', a.o + name, a.cpu, a.which))

        ## Concatenate ITS1 and ITS2 sequences if 'all' is selected
        if a.which == 'all':
            with open(a.o + name + '.ITS1.fasta', 'r') as handle1, open(a.o + name + '.ITS2.fasta', 'r') as handle2:
                its1_seqs = list(SeqIO.parse(handle1, 'fasta'))
                its2_seqs = list(SeqIO.parse(handle2, 'fasta'))
                
                concatenated_seqs = []
                for seq1, seq2 in zip(its1_seqs, its2_seqs):
                    concatenated_seq = SeqRecord(
                        seq1.seq + seq2.seq,
                        id = a.name,
                        name=seq1.name + "_concat_" + seq2.name,
                        description=seq1.description + " | " + seq2.description + " length: " + str(len(seq1.seq) + len(seq2.seq))
                    )
                    concatenated_seqs.append(concatenated_seq)
                
                with open(a.o + name + '.all.fasta', "w") as output_handle:
                    SeqIO.write(concatenated_seqs, output_handle, 'fasta')

        ## Remove Duplicates from concatenated or single ITSx output
        output_file = a.o + name + '.all.fasta' if a.which == 'all' else a.o + name + '.' + a.which + '.fasta'
        with open(output_file, 'r') as handle:
            its_seqs = list(SeqIO.parse(handle, 'fasta'))
        seqs = [str(s.seq) for s in its_seqs]
        diff = list(set(seqs))
        unique = []
        for n in range(len(diff)):
            if len(diff) == 1:
                Name = name
            else:
                Name = name + '_#' + str(n + 1)
            unique.append(SeqRecord(Seq(diff[n]), id=Name, name=Name, description=str(seqs.count(diff[n])) + ' copies of this sequence found in ' + a.i.split('/')[-1]))

        with open(a.o + name + '.' + a.which + '_filtered.fasta', "w") as output_handle:
            SeqIO.write(unique, output_handle, 'fasta')

        ## End message

        if len(diff) == 1:
            print('\nDONE. A single ' + a.which + ' sequence was found in ' + str(seqs.count(diff[0])) + ' copies, in file ' + a.i.split('/')[-1] + '.')
        elif len(diff) > 1:
            print('\nDONE. ' + str(len(diff)) + ' different ' + a.which + ' sequences were found in file ' + a.i.split('/')[-1] + '. See output files for more info.')
        elif len(diff) == 0:
            print('\nDONE. No ' + a.which + ' sequence found in file ' + a.i.split('/')[-1])

    else:
        print('\nDONE. No ' + a.which + ' sequence found in file ' + a.i.split('/')[-1] + '.')
        terminal('touch ' + a.o + name + '.' + a.which + '_filtered.fasta')