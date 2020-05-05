from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
import argparse
import math
import os
import time

sequence_of_interest = Seq("NDVTSLISTTYPYTGPPPMSHGSSTKYTLETIKRTYDYSRTSVEKTSKVFNIPRRKFCNCLEDKDELVKP", generic_protein)

def regular_search(seq, file_name):
    i = 0
    total_count = 0
    for seq_record in SeqIO.parse(file_name, "fasta"):
        total_count += 1
        index = seq_record.seq.find(seq)
        if (index != -1):
            i += 1
            print(f"Found Something!")
            print(f"ID: {seq_record.id}")
    print(total_count)

def split_file(file_name):
    num_lines = sum(1 for line in open(file_name))
    lines_per_file = math.ceil(num_lines / 10) #split into 10 files
    smallfile = None
    keep_going = False
    file_number = 0
    with open(file_name) as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0 or keep_going == True:
                if line.startswith(">"): #make sure to split on a new entry
                    if smallfile:
                        smallfile.close()
                    file_number += 1
                    small_filename = 'small_file_{}_{}'.format(file_name,file_number)
                    smallfile = open(small_filename, "w")
                    keep_going = False
                else:
                    keep_going = True
            smallfile.write(line)
        if smallfile:
            smallfile.close()

def cleanup_split_files():
    cwd = os.getcwd()
    files_in_directory = os.listdir(cwd)
    filtered_files = [file for file in files_in_directory if file.startswith("small_file")]
    for file in filtered_files:
        path_to_file = os.path.join(os.getcwd(), file)
        os.remove(path_to_file)

def download_file(file_name):
    #TODO write logic to download a file
    print("working on it")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sequence',help="Sequence of proteins ie. NDVTSL", type = str)
    parser.add_argument('-f', '--file', help="Specify file for which to search. It should be in the container or be a downloadable link", type = str, default = "influenza_xsmall.faa", required = False)
    parser.add_argument('-c', '--cut', help="Should we split the file into 10 files", dest = "cut", action = "store_true")
    parser.add_argument('-i','--iteration', help="If a cut was specified, specify the iteration of the file you want searched", type = int, required = False, default = 1)
    arguments = parser.parse_args()
    #split_file(arguments.file)
    #cleanup_split_files()
    if (arguments.file.startswith("http")):
        #assume we need to download it
        download_file(file_name = arguments.file)
        file_name = arguments.file.split("/")[-1]

    file_name = arguments.file
    if arguments.cut:
        print("nah")
        split_file(arguments.file)
        new_file_name = 'small_file_{}_{}'.format(file_name,arguments.iteration)

    start_time = time.time() #Timing execution
    regular_search(seq = arguments.sequence, file_name = file_name)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()