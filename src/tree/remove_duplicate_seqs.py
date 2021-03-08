#!/usr/bin/env python

import sys
import textwrap
from Bio import SeqIO

def seqs_parser(fasfile):
    rec_dic = {}
    for seq_record in SeqIO.parse(fasfile, "fasta"):
        if str(seq_record.seq) in rec_dic.keys(): continue
        else:
            rec_dic[str(seq_record.seq)] = str(seq_record.description)
    return rec_dic


def main():
    fasfile=sys.argv[1]
    ref = sys.argv[2]
    rec_dic = seqs_parser(fasfile)
    for key, value in rec_dic.items():
        if ref in value:
            value = value.split(' ')[0]
        print('>' + value)
        print(key)


if __name__ == '__main__':
    main()