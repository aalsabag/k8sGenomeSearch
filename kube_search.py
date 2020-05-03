from search import regular_search
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

import os
from kubernetes import client, config, utils
from kubernetes.client.rest import ApiException



sequence_of_interest = Seq("NDVTSLISTTYPYTGPPPMSHGSSTKYTLETIKRTYDYSRTSVEKTSKVFNIPRRKFCNCLEDKDELVKP", generic_protein)

def main():
    #create jobs
    #run jobs
    regular_search(sequence_of_interest)

if __name__ == '__main__':
    main()