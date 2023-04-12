import getopt
import sys
import os 
import pandas as pd
import numpy as np
from Bio import pairwise2
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# 1) gc
def gc(sequence):
    sequence = sequence.upper()
    return (sequence.count('G') + sequence.count('C')) / len(sequence)

#2) transcribe
def transcribe(sequence):
    return sequence.replace('T', 'U')

#3) reverse complement
def reverse_complement(sequence):
    sequence = sequence.upper()
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(sequence)])

#4) calc_nbases
def calc_nbases(sequence):
    sequence = sequence.upper()
    return sequence.count('N')

#5) is_valid
def is_valid_sequence(sequence,type):
    types = ['DNA','RNA',"protein"]
    type_sets = {'DNA': set('ATGCN'), 'RNA': set('AUGCN'), 'protein': set('ACDEFGHIKLMNPQRSTVWY')}
    if type not in types:
        raise ValueError("Invalid sequence type")
    return set(sequence.upper()).issubset(type_sets[type])

#6) filter nbases
def filter_nbases(sequence):
    return sequence.upper().replace('N', '')


def write_or_print(alignments , flag , output_file=None):
    if output_file:
        with open(output_file, 'w') as f:
            for alignment in alignments:
                score = alignment[2]
                seq1 = alignment[0]
                seq2 = alignment[1]
                f.write(f"Score -> {score}\n")
                f.write(f"{seq1}\n")
                f.write(f"{seq2}\n")
    else :
        for alignment in alignments:
            score = alignment[2]
            seq1 = alignment[0]
            seq2 = alignment[1]
            print(f"Score -> {score}")
            print(f"{seq1}")
            print(f"{seq2}")



# problem 7
def seq_alignment(seq1, seq2, output_file=None):
  # apply global alignment
  alignments = pairwise2.align.globalxx(seq1, seq2)

  # check if output_file != None
  if output_file:
    with open(output_file, 'w') as f:
        write_or_print(alignments, 1, 'output.txt')
  else:
        write_or_print(alignments, 0)
# seq_alignment("ATTGGCC", "AATGGCC" , 'output.txt')

# problem 8
def seq_alignment_files(file1, file2, output_file=None):
  # Read sequences from fasta files
  seq1 = SeqIO.read(file1, 'fasta').seq
  seq2 = SeqIO.read(file2, 'fasta').seq

  # apply global alignment
  alignments = pairwise2.align.globalxx(seq1, seq2)

  # check if output_file != None
  if output_file:
    with open(output_file, 'w') as f:
        write_or_print(alignments, 1, output_file)
  else:
        write_or_print(alignments, 0)
# seq_alignment_files(r'F:\sana 4\biopython\test\file1.fasta', r'F:\sana 4\biopython\test\file2.fasta')

# problem 9
def online_alignment(seq, output_file=None):
  # Perform BLAST search
  result_handle = NCBIWWW.qblast("blastn", "nr", seq)

  # Parse BLAST results
#   blast_record = Record.read(result_handle)
  blast_record = NCBIXML.read(result_handle)
  # Print or write results
  if output_file:
    with open(output_file, 'w') as f:
      for alignment in blast_record.alignments: 
        for hsp in alignment.hsps: 
            f.write('****Alignment****')
            f.write("\n")
            title = 'sequence:' + alignment.title
            f.write(title) 
            f.write(f'length:{alignment.length}') 
            f.write(f'e value: {hsp.expect}') 
            f.write(str(hsp.query)) 
            f.write(str(hsp.match)) 
            f.write(str(hsp.sbjct))
            f.write(str(hsp.score))
            f.write(str(hsp.bits))
            f.write(str(hsp.num_alignments))
            f.write(str(hsp.identities))
            f.write(str(hsp.positives))
            f.write(str(hsp.gaps))
            f.write(str(hsp.strand))
            f.write(str(hsp.frame))
            f.write(str(hsp.query_start))
            f.write(str(hsp.sbjct_start))
  else:
    for alignment in blast_record.alignments: 
        for hsp in alignment.hsps: 
            print('****Alignment****') 
            print('sequence:', alignment.title) 
            print('length:', alignment.length) 
            print('e value:', hsp.expect) 
            print(hsp.query) 
            print(hsp.match) 
            print(hsp.sbjct)
            print(hsp.score)
            print(hsp.bits)
            print(hsp.num_alignments)
            print(hsp.identities)
            print(hsp.positives)
            print(hsp.gaps)
            print(hsp.strand)
            print(hsp.frame)
            print(hsp.query_start)
            print(hsp.sbjct_start)

# search_alignments('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAACGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGTGACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCCCGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCCCAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAACGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTGAATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCAGGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCGGCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCGGCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTGGCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCCTTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGGGGCACCCGCTGAGTTTACGC','output.txt')

# problem 10
def merge_fasta(*fasta_files, output_file=None):
    total_seq = []
    for fasta_file in fasta_files:
        all_records = SeqIO.parse(fasta_file, 'fasta')
        for record in all_records:
            total_seq.append(record)
        print (total_seq)
    if output_file:
        SeqIO.write(total_seq, output_file, 'fasta')
    else:
        for sequence in total_seq:
            print(f">{sequence.id}\n{sequence.seq}")

# merge_fasta(r'F:\sana 4\biopython\test\file1.fasta', r'F:\sana 4\biopython\test\file2.fasta')

# problem 11
def convert_to_fasta(gb_file):
    output_file = "output.fasta"
    total_seq = []
    for record in SeqIO.parse(gb_file, "genbank"):
        total_seq.append(record)
    SeqIO.write(total_seq, output_file, 'fasta')
    
# convert_to_fasta(r"F:\sana 4\biopython\test\ls_orchid.gbk")

commands_list = ["gc" ,
                 "transcribe",
                 "reverse_complement",
                 "calc_nbases",
                 "is_valid",
                 "filter_nbases",
                 "seq_alignment",
                 "seq_alignment_files",
                 "online_alignment",
                 "merge_fasta",
                 "convert_to_fasta",
                 "help"]

code_name = sys.argv[0]
command = sys.argv[1]

if command in commands_list:
    # 1. gc
    if command == "gc":
        if len(sys.argv) == 3:
            x = gc(sys.argv[2])
            print(x)
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 2.transcribe
    elif command == "transcribe":
        if len(sys.argv) == 3:
            x=transcribe(sys.argv[2])
            print(x)
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 3. reverse_complement
    elif command == "reverse_complement":
        if len(sys.argv) == 3:
            x=reverse_complement(sys.argv[2])
            print(x)
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 4. calc_nbases
    elif command == "calc_nbases":
        if len(sys.argv) == 3:
            x=calc_nbases(sys.argv[2])
            print(x)
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 5. is_valid
    elif command == "is_valid":
        if len(sys.argv) == 3:
            x=is_valid_sequence(sys.argv[2])
            if x:
                print("True")
            else:
                print("False")
        
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 6. filter_nbases
    elif command == "filter_nbases":
        if len(sys.argv) == 3:
            x=filter_nbases(sys.argv[2])
            print(x)
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 7. seq_alignment
    elif command == "seq_alignment":
        opts , args = getopt.getopt(sys.argv[4:], "o:",['output ='])
        if len(sys.argv) == 6:
            for opt, arg in opts:
                if opt in ("-o", "--output"):
                    seq_alignment(sys.argv[2], sys.argv[3], arg)
        elif len(sys.argv) == 4:
            seq_alignment(sys.argv[2], sys.argv[3])
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 8. seq_alignment_files
    elif command == "seq_alignment_files":
        opts , args = getopt.getopt(sys.argv[4:], "o:",["output"])
        print(opts)
        print(args)
        print(sys.argv)
        if len(sys.argv) == 6:
            for opt, arg in opts:
                if opt in ("-o", "--output"):
                    seq_alignment_files(sys.argv[2], sys.argv[3], arg)
        elif len(sys.argv) == 4:
            seq_alignment_files(sys.argv[2], sys.argv[3])
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 9. online_alignment
    elif command == "online_alignment":
        opts , args = getopt.getopt(sys.argv[3:], "o:",["output"])
        if len(sys.argv) == 5:
            for opt, arg in opts:
                if opt in ("-o", "--output"):
                    online_alignment(sys.argv[2], arg)
        elif len(sys.argv) == 3:
            online_alignment(sys.argv[2])
    # 10. merge_fasta
    elif command == "merge_fasta":
        opts , args = getopt.getopt(sys.argv[-2:], "o:",["output"])
        print("hello")
        print(opts)
        print(args)
        print(sys.argv)
        if len(opts) > 0:
            for opt, arg in opts:
                if opt in ("-o", "--output"):
                    merge_fasta(*sys.argv[2:-2], output_file=arg)
        elif len(opts) == 0:
            merge_fasta(*sys.argv[2:])
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    # 11. convert_to_fasta
    elif command == "convert_to_fasta":
        if len(sys.argv) == 3:
            convert_to_fasta(sys.argv[2])
        else:
            print("Invalid command. Please enter a valid command. Type help for list of commands.")
            exit()
    
    else:
        print("Invalid command. Please enter a valid command. Type help for list of commands.")
        exit()