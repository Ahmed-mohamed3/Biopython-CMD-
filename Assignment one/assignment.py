from Bio import pairwise2
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


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
        write_or_print(alignments, 1, 'output.txt')
  else:
        write_or_print(alignments, 0)
# seq_alignment_files(r'F:\sana 4\biopython\test\file1.fasta', r'F:\sana 4\biopython\test\file2.fasta')

# problem 9
def search_alignments(seq, output_file=None):
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
            f.write('sequence:', alignment.title) 
            f.write('length:', alignment.length) 
            f.write('e value:', hsp.expect) 
            f.write(hsp.query) 
            f.write(hsp.match) 
            f.write(hsp.sbjct)
            f.write(hsp.score)
            f.write(hsp.bits)
            f.write(hsp.num_alignments)
            f.write(hsp.identities)
            f.write(hsp.positives)
            f.write(hsp.gaps)
            f.write(hsp.strand)
            f.write(hsp.frame)
            f.write(hsp.query_start)
            f.write(hsp.sbjct_start)
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
def gconvert_to_fasta(gb_file, output_file=None):
    total_seq = []
    if output_file:
        for record in SeqIO.parse(gb_file, "genbank"):
            total_seq.append(record)
        SeqIO.write(total_seq, output_file, 'fasta')
    else:
        all_records = list(SeqIO.parse(gb_file, "genbank"))
        for sequence in all_records:
            print(f">{sequence.id}\n{sequence.seq}")
# convert_to_fasta(r"F:\sana 4\biopython\test\ls_orchid.gbk",'output.fasta')