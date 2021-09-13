#!/usr/bin/env python

__author__ = 'atamuri'

"""
./prune_alignment.py <alignment_file> <file_of_taxa_to_remove>
"""

from Bio import SeqIO, AlignIO
import sys
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filein', type = argparse.FileType('r'), default = sys.stdin, nargs='?')
    parser.add_argument('--keep', type = argparse.FileType('r'), nargs=1)
    parser.add_argument('--remove', type = argparse.FileType('r'), nargs=1)
    args = parser.parse_args()
    args = vars(args)

    if args['remove'] is None and args['keep'] is None:
        print "Must supply a file of taxa to --remove <file> or --keep <file>"
        print parser.print_help()
        sys.exit(1)
    elif args['remove'] is not None and args['keep'] is not None:
        print "Must pick --remove <file> OR --keep <file>"
        print parser.print_help()
        sys.exit(1)

    seqs = list(SeqIO.parse(args['filein'], format="fasta"))

    if args['remove'] is not None:
        remove = args['remove'][0].read().splitlines()
        pruned = [s for s in seqs if s.id not in remove]
    elif args['keep'] is not None:
        keep = args['keep'][0].read().splitlines()
        pruned = [s for s in seqs if s.id in keep]

    alignment = AlignIO.MultipleSeqAlignment(pruned)

    handle = sys.stdout
    AlignIO.write(alignment, handle, "fasta")
    sys.stdout.flush()

