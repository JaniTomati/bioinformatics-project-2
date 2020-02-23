#!/usr/bin/python3
# -*- coding: utf-8 -*-

import time
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Example call
# python3 global_alignment_linear_gapcost.py --C scoring_matrix.txt --P Y Seq1.fasta Seq2.fasta

# Default values
cost = np.array([
    # A  C  G  T
    [10, 2, 5, 2],
    [2, 10, 2, 5],
    [5, 2, 10, 2],
    [2, 5, 2, 10]])
gap_cost = -5
look_up = {"A": 0, "C": 1, "G": 2, "T":3}


# Setting parameters
def read_in_scoring_matrix(file):
    """ Read in the alphabet, scoring matrix and gap cost """
    alphabet = []
    scoring_matrix = []

    with open(file, "r") as f:
        lines = f.readlines()

    firstLine = True
    for line in lines:
        line_arr = line.rstrip().split()
        if firstLine == True:
            gapcost = int(line_arr[0])
            firstLine = False
        else:
            row = []
            for char in line_arr:
                if char == line_arr[0]: # add char to alphabet
                    alphabet.append(char)
                else:                   # add cost to scoring matrix
                    row.append(int(char))
            scoring_matrix.append(row)

    return gapcost, np.array(scoring_matrix), alphabet


# Optimal alignment
def optimal_alignment(seq1, seq2):
    m, n = len(seq1)+1, len(seq2)+1
    align_matrix = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            if j == 0 and i == 0:
                align_matrix[i, j] = 0
            elif j == 0:
                align_matrix[i, j] = align_matrix[i-1, j] + gap_cost
            elif i == 0:
                align_matrix[i, j] = align_matrix[i, j-1] + gap_cost

            else:
                val = min(
                    align_matrix[i-1, j-1] + cost[look_up[seq2[i-1]], look_up[seq1[j-1]]],
                    align_matrix[i-1, j] + gap_cost,
                    align_matrix[i, j-1] + gap_cost
                )
                align_matrix[i, j] = val

    return align_matrix


# Traceback algorithm
def traceback(opt, seq1, seq2):
    s1 = []
    s2 = []
    i, j = len(seq2), len(seq1)
    while i > 0 and j > 0:
        if opt[i, j] == opt[i-1, j-1]+cost[look_up[seq1[j-1]], look_up[seq2[i-1]]]:
            s2.append(seq2[i-1])
            s1.append(seq1[j-1])
            i -= 1
            j -= 1

        elif opt[i-1, j] + gap_cost == opt[i,j]:
            s1.append("-")
            s2.append(seq2[i-1])
            i -= 1

        elif opt[i, j-1] + gap_cost == opt[i,j]:
            s2.append("-")
            s1.append(seq1[j-1])
            j -= 1

        elif i == 0:
            while j >= 0:
                s2.append("-")
                j -= 1
        elif j == 0:
            while i >= 0:
                s1.append("-")
                i -= 1

    s1 = "".join(s1[::-1])
    s2 = "".join(s2[::-1])

    return s1, s2


def measure_times():
    """ Use this method to measure times of the algorithms using sequences of different length """
    sequences = []
    sequence_length = [5, 10, 50, 100, 500, 1000, 2500, 5000, 7500, 10000] # edit sequence

    for seq_record in SeqIO.parse("sequences/rand_sequences.fasta", "fasta"):
        sequences.append(seq_record.seq.upper())

    times_opt = []
    times_back = []
    i = 0
    while i < len(sequences):
        # measure optimal alignment algorithm
        start_opt = time.time()
        opt = optimal_alignment(sequences[i], sequences[i + 1])
        end_opt = time.time()
        times_opt.append(end_opt - start_opt)

        # measure trace back algorithm
        start_back = time.time()
        aligned = traceback(opt, sequences[i], sequences[i + 1])
        end_back = time.time()
        times_back.append(end_back - start_back)
        i += 2

    return times_opt, times_back


def test_alignment_algorithm():
    """ Test the alignment algorithm by using the sequences from the project2_examples.txt """
    sequences = []

    for seq_record in SeqIO.parse("sequences/test_sequences.fasta", "fasta"):
        sequences.append((seq_record.seq).upper())

    i = 0
    while i < len(sequences):
        opt = optimal_alignment(sequences[i], sequences[i + 1])
        print("\nScore of the optimal global optimal_alignment: ", opt[len(sequences[i + 1]), len(sequences[i])])

        aligned = traceback(opt, sequences[i], sequences[i + 1])
        print(aligned[0], "\n" + aligned[1], "\n")
        i += 2


def main():
    global gap_cost, cost
    # Arg parse
    parser = argparse.ArgumentParser()
    parser.add_argument("Seq1", help = "Fasta file path for sequence 1")
    parser.add_argument("Seq2", help = "Fasta file path for sequence 2")
    parser.add_argument("--P", help = "Output strings of aligned sequences: Write Y")
    parser.add_argument("--C", help = "txt file with Phylip-like format")
    args = parser.parse_args()
    seq1 = args.Seq1
    seq2 = args.Seq2
    print_statement = args.P
    file = args.C

    if print_statement is not None:
        print_statement = print_statement.lower() # lower case letters

    if file is not None:
        gap_cost, cost, alphabet = read_in_scoring_matrix(file)

    ### Running the Code ###
    seq_list1 = []
    for seq_record in SeqIO.parse(seq1, "fasta"):
        seq_list1.append((seq_record.seq))

    seq_list2 = []
    for seq_record in SeqIO.parse(seq2, "fasta"):
        seq_list2.append((seq_record.seq))

    seq1 = seq_list1[0]
    seq2 = seq_list2[0]

    start_opt = time.time()
    opt = optimal_alignment(seq1, seq2)
    end_opt = time.time()
    print("\nScore of the optimal global optimal_alignment: ", opt[len(seq2), len(seq1)])
    print("Calculating the optimal alignment took", end_opt - start_opt, "seconds.\n")

    if print_statement == "y":
        start_back = time.time()
        aligned = traceback(opt, seq1, seq2)
        end_back = time.time()

        #print("",aligned[0], "\n", aligned[1])
        print(aligned[0])
        print(aligned[1], "\n")
        print("Calculating the traceback took", end_back - start_back, "seconds.\n")

    ### Writing the aligned sequences into a fasta file ###
    record = SeqRecord(Seq(aligned[0]),
                       id="Seq1")
    record1 = SeqRecord(Seq(aligned[1]),
                       id="Seq2")

    result_file = "result_linear.fasta"
    SeqIO.write([record, record1], result_file, "fasta")
    print("Wrote results in fasta format to " +  result_file + ".")

    # test alignment algorithm using examples
    test_alignment_algorithm()

    # function for taking the performance measures
    # times_opt, times_back = measure_times()
    # print("Optimal alignment:", times_opt)
    # print("Traceback", times_back)


if __name__ == '__main__':
    main()
