from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
import argparse

sequence_of_interest = Seq("NDVTSLISTTYPYTGPPPMSHGSSTKYTLETIKRTYDYSRTSVEKTSKVFNIPRRKFCNCLEDKDELVKP", generic_protein)

def regular_search(seq):
    i = 0
    total_count = 0
    for seq_record in SeqIO.parse("influenza_xsmall.faa", "fasta"):
        total_count += 1
        index = seq_record.seq.find(seq)
        if (index != -1):
            i += 1
            print(f"Found Something. Item #{i}")
            print(f"ID: {seq_record.id}")
            print(f"Index in the sequence: {index}")
    print(total_count)
            

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sequence',help="Sequence of proteins ie. NDVTSL", type = str)
    arguments = parser.parse_args()
    regular_search(seq = arguments.sequence)

if __name__ == '__main__':
    main()