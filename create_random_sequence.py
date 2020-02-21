#!/usr/bin/python3
# -*- coding: utf-8 -*-

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def create_random_sequence(n):
    """ Create a random DNA sequence of length n """
    alphabet = ["A", "C", "G", "T"]
    sequence = ""
    for i in range(n + 1):
        sequence += random.choice(alphabet)

    return sequence


def main():
    n = [5, 10, 100, 1000, 10000, 100000]
    sequences = []

    count = 0
    index = "Seq"
    for sequence_length in n:
        sequence = create_random_sequence(sequence_length)
        sequences.append(SeqRecord(Seq(sequence, generic_dna), id=(index + str(count)), description=""))
        count += 1

        sequence = create_random_sequence(sequence_length)
        sequences.append(SeqRecord(Seq(sequence, generic_dna), id=(index + str(count)), description=""))
        count += 1

    with open("sequences/rand_sequences.fasta", "w") as out_file:
        SeqIO.write(sequences, out_file, "fasta")



if __name__ == '__main__':
    main()
