#!/usr/bin/python3
# -*- coding: utf-8 -*-

import time
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Default values
cost = np.array([
    # A  C  G  T
    [10, 2, 5, 2],
    [2, 10, 2, 5],
    [5, 2, 10, 2],
    [2, 5, 2, 10]])
gapcost = (5, 5) # alpha, beta
look_up = {"A": 0, "C": 1, "G": 2, "T":3}


def read_in_scoring_matrix(file):
    """ Read in the alphabet, scoring matrix and gap cost """
    gapcost = None
    scoring_matrix = []
    alphabet = []

    with open(file, "r") as f:
        lines = f.readlines()

    firstLine = True
    for line in lines:
        line_arr = line.rstrip().split()
        if firstLine == True:
            if len(line_arr) != 1:
                gapcost = (int(line_arr[0]), int(line_arr[1]))
            else:
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

    return gapcost, alphabet, np.array(scoring_matrix)


# Affine Alignment
def affine_alignment(seq1, seq2, gap_start=5, gap_extend=5):
    n=len(seq1) + 1
    m=len(seq2) + 1
    S, D, I = np.zeros((n,m)), np.zeros((n,m)), np.zeros((n,m))

    for j in range(0,m):
        for i in range(0,n):
            if i==0 and j==0:
                S[i,j] = 0
                I[i,j] = gap_start
                D[i,j] = gap_start
            elif i==0:
                I[i,j] = I[i,j-1] + gap_extend
                S[i,j] = I[i,j]
                D[i,j] = S[i,j]
            elif j==0:
                D[i,j] = D[i-1,j] + gap_extend
                S[i,j] = D[i,j]
                I[i,j] = S[i,j]
            else:
                I[i,j]=min(I[i,j-1] + gap_extend,
                           S[i,j-1] + gap_start + gap_extend)
                D[i,j]=min(D[i-1,j] + gap_extend,
                           S[i-1,j] + gap_start+gap_extend)
                S[i,j]=min(I[i,j],
                           D[i,j],
                           S[i-1,j-1] + cost[look_up[seq1[i-1]], look_up[seq2[j-1]]])

    return S


def penalty(k, gap_extend=5, gap_start = 5):
        return gap_start+k*gap_extend


def traceback(matrix, seq1, seq2):
    i = len(seq1)
    j = len(seq2)
    align1 = ""
    align2 = ""
    while i>0 and j>0:
        if matrix[i,j] == matrix[i-1,j-1] + cost[look_up[seq1[i-1]], look_up[seq2[j-1]]]:
            align1 = str(seq1[i-1]) + align1
            align2 = str(seq2[j-1]) + align2
            i -= 1
            j -= 1
        else:
            k=1
            while True:
                if matrix[i,j] == matrix[i,j-k] + penalty(k):
                    align1 = "-" * k + align1
                    align2 = str(seq2[j-k:j]) + align2
                    j -= k
                    break
                elif matrix[i,j] == matrix[i-k,j] + penalty(k):
                    align1 = str(seq1[i-k:i])+align1
                    align2 = "-" * k + align2
                    i -= k
                    break
                else:
                    k += 1

    return align1, align2


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
        print_statement = print_statement.lower() # so it is allways a small letter

    if file is not None:
        gap_cost, alphabet, cost = read_in_scoring_matrix(file)

    ### Running the Code ###
    seq_list1 = []
    for seq_record in SeqIO.parse(seq1, "fasta"):
        seq_list1.append((seq_record.seq))

    seq_list2 = []
    for seq_record in SeqIO.parse(seq2, "fasta"):
        seq_list2.append((seq_record.seq))

    seq1 = seq_list1[0]
    seq2 = seq_list2[0]

    start_aff = time.time()
    aff = affine_alignment(seq1, seq2, gapcost[0], gapcost[1])
    end_aff = time.time()
    print("\n", "Score of the optimal global optimal_alignment: ", aff[len(seq1), len(seq2)], "\n")
    print("Calculating the optimal alignment took", end_aff - start_aff, "seconds.")

    if print_statement == "y":
        start_back = time.time()
        aligned = traceback(aff, seq1, seq2)
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

    result_file = "result_affine.fasta"
    SeqIO.write([record, record1], result_file, "fasta")
    print("Wrote results in fasta format to " +  result_file + ".")


if __name__ == '__main__':
    main()
