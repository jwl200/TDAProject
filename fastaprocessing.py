from __future__ import print_function
import sys, re
import os
import numpy as np
import subprocess, fileinput

'''Computes hamming distance matrix from input FASTA file'''

def Hamming_distance(string1, string2):
    '''Computes the Hamming distances between two sequences for a given alignment
    of genome sequences.'''
    if len(string1) != len(string2):
        raise ValueError("Sequences must be of the same length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(string1, string2) if ch1 != "N" and ch2 != "N")

def format_data(input_file, filetype):
    '''Formats input data for Hamming distance matrix computations'''
    if filetype != "fasta":
        with open(input_file, 'r') as f:
            lines = f.readlines()
        for i in xrange(len(lines)):
            lines[i] = lines[i][1:]
    else:
        out = []
        with open(input_file, 'r') as f:
            lines = f.readlines()[1::2]
        x = min([len(i) for i in lines])
        #print(x)            
    for line in lines:
        out.append(line[0:x])
    return out

def single_line_fasta(fasta_in, fasta_out): #moves each sequece to a single line
    with open(fasta_in) as f_input, open(fasta_out, 'w') as f_output:
        block = []
    
        for line in f_input:
            if line.startswith('>'):
                if block:
                    f_output.write(''.join(block) + '\n')
                    block = []
                f_output.write(line)
            else:
                block.append(line.strip())
    
        if block:
            f_output.write(''.join(block) + '\n')
            
def empty_matrix(lines):
    '''Generate an empty matrix to fill with Hamming distances'''
    return np.ndarray(shape = (len(lines), len(lines)), dtype = float)

def populate_matrix(matrix, lines):
    '''Populates matrix elements with pairwise Hamming distances'''
    for i in range(len(lines)):
        for j in range(len(lines)):
            matrix[i][j] = Hamming_distance(lines[i], lines[j])
    return matrix

