#!/usr/bin/env python

#Given a sequence, this program uses a Random Forest to predict the activity
#via a feature matrix, using a predefined master model created using a training set

# USAGE:
# Users can use the TUSCAN model in 3 ways:
# 1. Only supply a fasta file. The entire sequence will be analysed. Output will show: Chrom: N/A. Candidate site locations will be 0-based.
# 2. Supply a fasta file with known chrom/position information. The entire sequence will be analysed. Output will show supplied chromosome info, and position relative to supplied start/finish positions. 
# 3. Supply a fasta file, include chrom/position information for a region to extract. Include the -e flag with no arguments. The region of interest (start pos -> end pos) will be extracted and analysed. Output will reflect start position relative to supplied start/end position and chrom info.

# Users can also supply a thread number for multithreading using -t. If none supplied, default thread depending on available cores will be set. (-t [THREADS (int)])
# Users can also supply an output file name (recommended), otherwise output is written to TUSCAN_output.txt and overwritten upon each iteration (-o [FILENAME])
# Users must also specify whether regression or classification is required (-m [Regression/Classification])

# NOTE: Final output will be in arbitrary order due to multithreading. 
# Output is tab-separated and easy for users to manipulate using UNIX sort or other flavours. 
# Ouput is in bed6 format, with the unique ID being chrom_start_stop_targetSequence

import sys
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.externals import joblib
import numpy
import pybedtools
import argparse
from collections import OrderedDict, namedtuple
import os
import re
# from string import maketrans
from multiprocessing import Process, Queue, cpu_count

NucleotideLocation = namedtuple('NucleotideLocation', ['nucleotide', 'location'])

DinucleotideLocation = namedtuple('DinucleotideLocation', ['dinucleotide', 'location'])

REGRESSION_NUCLEOTIDES_OF_INTEREST = (
    NucleotideLocation(nucleotide='T', location=4),
    NucleotideLocation(nucleotide='C', location=7),
    NucleotideLocation(nucleotide='C', location=10),
    NucleotideLocation(nucleotide='T', location=17),
    NucleotideLocation(nucleotide='C', location=20),
    NucleotideLocation(nucleotide='T', location=20),
    NucleotideLocation(nucleotide='G', location=21),
    NucleotideLocation(nucleotide='T', location=21),
    NucleotideLocation(nucleotide='G', location=22),
    NucleotideLocation(nucleotide='T', location=22),
    NucleotideLocation(nucleotide='C', location=24),
    NucleotideLocation(nucleotide='G', location=24),
    NucleotideLocation(nucleotide='T', location=24)
)

CLASSIFICATION_NUCLEOTIDES_OF_INTEREST = (
    NucleotideLocation(nucleotide='G', location=5),
    NucleotideLocation(nucleotide='T', location=11),
    NucleotideLocation(nucleotide='C', location=12),
    NucleotideLocation(nucleotide='A', location=16),
    NucleotideLocation(nucleotide='T', location=17),
    NucleotideLocation(nucleotide='C', location=20),
    NucleotideLocation(nucleotide='T', location=20),
    NucleotideLocation(nucleotide='T', location=22),
    NucleotideLocation(nucleotide='T', location=23),
    NucleotideLocation(nucleotide='C', location=24),
    NucleotideLocation(nucleotide='G', location=24),
    NucleotideLocation(nucleotide='T', location=24),
)

REGRESSION_DINUCLEOTIDES_OF_INTEREST = (
    DinucleotideLocation(dinucleotide='AC', location=1),
    DinucleotideLocation(dinucleotide='AC', location=2),
    DinucleotideLocation(dinucleotide='CA', location=3),
    DinucleotideLocation(dinucleotide='TT', location=4),
    DinucleotideLocation(dinucleotide='GA', location=5),
    DinucleotideLocation(dinucleotide='CT', location=6),
    DinucleotideLocation(dinucleotide='AC', location=8),
    DinucleotideLocation(dinucleotide='CC', location=8),
    DinucleotideLocation(dinucleotide='GA', location=8),
    DinucleotideLocation(dinucleotide='TT', location=9),
    DinucleotideLocation(dinucleotide='AT', location=10),
    DinucleotideLocation(dinucleotide='CG', location=11),
    DinucleotideLocation(dinucleotide='GA', location=12),
    DinucleotideLocation(dinucleotide='CC', location=14),
    DinucleotideLocation(dinucleotide='GA', location=15),
    DinucleotideLocation(dinucleotide='CC', location=16),
    DinucleotideLocation(dinucleotide='GG', location=16),
    DinucleotideLocation(dinucleotide='TT', location=16),
    DinucleotideLocation(dinucleotide='CT', location=17),
    DinucleotideLocation(dinucleotide='AA', location=18),
    DinucleotideLocation(dinucleotide='GG', location=19),
    DinucleotideLocation(dinucleotide='AT', location=20),
    DinucleotideLocation(dinucleotide='CC', location=20),
    DinucleotideLocation(dinucleotide='CG', location=20),
    DinucleotideLocation(dinucleotide='CT', location=20),
    DinucleotideLocation(dinucleotide='GG', location=20),
    DinucleotideLocation(dinucleotide='TA', location=21),
    DinucleotideLocation(dinucleotide='TG', location=21),
    DinucleotideLocation(dinucleotide='CC', location=22),
    DinucleotideLocation(dinucleotide='GA', location=22),
    DinucleotideLocation(dinucleotide='TA', location=22),
    DinucleotideLocation(dinucleotide='CG', location=23),
    DinucleotideLocation(dinucleotide='GA', location=23),
    DinucleotideLocation(dinucleotide='GG', location=23),
    DinucleotideLocation(dinucleotide='TG', location=23),
    DinucleotideLocation(dinucleotide='GA', location=24),
    DinucleotideLocation(dinucleotide='GT', location=24),
    DinucleotideLocation(dinucleotide='TC', location=24),
)

CLASSIFICATION_DINUCLEOTIDES_OF_INTEREST = (
    DinucleotideLocation(dinucleotide='CG', location=11),
    DinucleotideLocation(dinucleotide='GA', location=15),
    DinucleotideLocation(dinucleotide='TT', location=16),
    DinucleotideLocation(dinucleotide='CC', location=20),
    DinucleotideLocation(dinucleotide='TA', location=22),
    DinucleotideLocation(dinucleotide='CG', location=23),
    DinucleotideLocation(dinucleotide='TC', location=23),
    DinucleotideLocation(dinucleotide='TG', location=23),
    DinucleotideLocation(dinucleotide='CC', location=24),
    DinucleotideLocation(dinucleotide='GA', location=24),
    DinucleotideLocation(dinucleotide='GC', location=24),
    DinucleotideLocation(dinucleotide='GT', location=24),
    DinucleotideLocation(dinucleotide='TC', location=24),
)

CLASSIFICATION_DINUCLEOTIDES = [
    'AA',
    'AC',
    'AG',
    'AT',
    'CA',
    'CC',
    'CG',
    'CT',
    'GA',
    'GC',
    'GG',
    'GT',
    'TA',
    'TC',
    'TG',
    'TT',
]

REGRESSION_DINUCLEOTIDES = [
    'CA',
    'CT',
    'GC',
    'TC',
    'TG',
    'TT',
]

#determines gc content of given sequence
def gc(seq, features, index):
    features[index] = round((seq.count('C') + seq.count('G'))/float(len(seq)) * 100, 2)

#counts appearance of dinucleotides in sequence
def di_content(seq, dinucleotides_to_count, features, start_index):
    for idx, dinucleotide in enumerate(dinucleotides_to_count):
        count = start = 0
        while True:
            start = seq.find(dinucleotide, start) + 1
            if start:
                count += 1
            else:
                features[start_index+idx] = count
                break

#checks if specific PAM is present in sequence
def pam(seq, features, index):
    if seq[24:28] == 'TGGT':
        features[index] = 1

#checks if given position-specific nucleotides are present in sequence
def nucleotide(seq, nucleotides_of_interest, features, start_index):
    for idx, nucleotide_loc in enumerate(nucleotides_of_interest):
        if seq[nucleotide_loc.location-1] == nucleotide_loc.nucleotide:
            features[start_index+idx] = 1

#checks if given position-specific dinucleotides are present in sequence
def dinucleotide(seq, dinucleotides_of_interest, features, start_index):
    #-1 is since a sequence of length N has N-1 dinucleotides
    for idx, dinucleotides_loc in enumerate(dinucleotides_of_interest):
        location = dinucleotides_loc.location
        if seq[location-1:location+1] == dinucleotides_loc.dinucleotide:
            features[start_index+idx] = 1

#generates a feature vector from a given 30 nucleotide sequence
def get_features(seq, is_regression):
    if is_regression:
        features = [0] * 63
        gc(seq, features, 0)
        features[1] = seq.count('A')
        features[2] = seq.count('C')
        features[3] = seq.count('G')
        features[4] = seq.count('T')
        di_content(seq, REGRESSION_DINUCLEOTIDES, features, 5)
        nucleotide(seq, REGRESSION_NUCLEOTIDES_OF_INTEREST, features, 11)
        dinucleotide(seq, REGRESSION_DINUCLEOTIDES_OF_INTEREST, features, 24)
        pam(seq, features, 62)
    else:
        features = [0] * 46
        gc(seq, features, 0)
        features[1] = seq.count('A')
        features[2] = seq.count('C')
        features[3] = seq.count('G')
        features[4] = seq.count('T')
        di_content(seq, CLASSIFICATION_DINUCLEOTIDES, features, 5)
        nucleotide(seq, CLASSIFICATION_NUCLEOTIDES_OF_INTEREST, features, 21)
        dinucleotide(seq, CLASSIFICATION_DINUCLEOTIDES_OF_INTEREST, features, 33)
    return features

def output_sequences(sequences, feature_lists, output_queue):
    feature_array = numpy.array(feature_lists)
    scores = rf.predict(feature_array)
    # output_queue.put('\n'.join(seq + '\t' + str(scores[i]) for i, seq in enumerate(sequences)))
    output_queue.put('\n'.join(seq[0] + '\t' + seq[1] + '\t' + seq[2] + '\t' + seq[3] + '\t' + str(scores[i]) + '\t' + seq[4] for i, seq in enumerate(sequences)))

def score_sequences(matches_queue, output_queue, is_reverse):
    strand = "-" if is_reverse else "+"
    sequences = []
    feature_lists = []
    while True:
        match_start, sequence = matches_queue.get()
        if match_start == "EMPTY":
            if sequences:
                output_sequences(sequences, feature_lists, output_queue)
            output_queue.put("DONE")
            break
        if is_reverse:
            sequence_end_pos = end - match_start - 3 - 1
            sequence_start_pos = sequence_end_pos - 23 + 1 
        else:
            sequence_start_pos = start + match_start + 3 + 2
            sequence_end_pos = sequence_start_pos + 23 - 1
        feature_list = get_features(sequence, is_regression)
        # sequences.append('{!s:5} {!s:10} {!s:10} {!s:8} {!s:31}'.format(chrom, str(sequence_start_pos), str(sequence_end_pos), strand, sequence[4:-3]))
        uniq_id = str(str(chrom)+"_"+str(sequence_start_pos)+"_"+str(sequence_end_pos)+"_"+str(sequence[4:-3]))
        sequences.append([chrom, str(sequence_start_pos), str(sequence_end_pos), uniq_id, strand])
        feature_lists.append(feature_list)
        if len(sequences) >= 10000:
            output_sequences(sequences, feature_lists, output_queue)
            sequences = []
            feature_lists = []

def fill_queue(matches, matches_queue):
    while True:
        try:
            match = next(matches)
        except StopIteration:
            for x in range(num_threads):
                matches_queue.put(("EMPTY", "EMPTY"), block=True)
            break
        matches_queue.put((match.start(), match.group(1)), block=True)

#collect directory
dir = os.path.dirname(os.path.realpath(__file__))

try:
    num_cores = cpu_count()
except NotImplementedError:
    num_cores = 4

#Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', required=False, help='Genome')
parser.add_argument('-o', required=False, help='Output file')
parser.add_argument('-s', required=False, help='Start Position')
parser.add_argument('-f', required=False, help='End Position')
parser.add_argument('-c', required=False, help='Chromosome')
parser.add_argument('-m', required=True, help='Type of Model: (Regression/Classification)')
parser.add_argument('-t', required=False, default=num_cores, help='Number of threads (default: 4)')
parser.add_argument('-e', required=False, action='store_true', help='If you want to excise a region inside a supplied genome, include this flag with no arguments, exclude it to analyse the full genome')

args = parser.parse_args()

sequence = ""
gaveGenome = False
is_regression = False
if args.m == "Regression":
    is_regression = True
elif args.m != "Classification":
    print("Please specify Regression or Classification")
    sys.exit()

chrom = args.c if args.c else "N/A"
genome = args.g
output_file = args.o if args.o else "TUSCAN_output.txt"
num_threads = args.t
extract = args.e
#if a chromosome, genome, start and stop positions are given
if (extract):
    print("Extracting region")
    #correct region info supplied
    gaveGenome = True
    start = int(args.s)
    end = int(args.f)
    #turn region information into string and write to BED file
    region = str(chrom + "\t" + str(start) + "\t" + str(end))
    f = open("customRegion.bed", "w")
    f.write(region)
    f.close()
    #extract region of interest from genome
    c = pybedtools.BedTool("customRegion.bed")
    d = c.sequence(fi = genome)
    with open(d.seqfn, 'r') as g:
        for line in g:
            line = line.rstrip()
            if not line.startswith(">"):
                sequence = str(line).upper()
    os.remove("customRegion.bed")
#else if a FASTA sequence has been nominated
elif (genome):
    #extract from file:
    print("Analysing entire supplied sequence.\nIf you wish to analyse a sub-region, supply start, end and chromosome flags and include the -e flag")
    sequence = ""
    with open(genome, 'r') as g:
        for line in g:
            line = line.rstrip()
            if not line.startswith(">"):
                sequence += str(line).upper()
    start = int(args.s) - 1 if args.s else 0
    end = int(args.f) - 1 if args.f else len(sequence)
else:
    print("Please supply a sequence to be analysed")
    sys.exit()

if not sequence:
    print("If using the extract -e flag, please supply ALL OF THESE and make sure they are correct: start, end AND chromosome flags")
    sys.exit()

LAYOUT = '{!s:5} {!s:10} {!s:10} {!s:8} {!s:34} {!s:15}'
header = LAYOUT.format('Chrom', 'Start', 'End', 'ID', 'TUSCANScore', 'Strand')

#Find and store all sequences + location information
input_base = "ATCG"
output_base = "TAGC"
complement = str.maketrans(input_base, output_base)
reverse_sequence = sequence.translate(complement)[::-1]

matches = re.finditer(r'(?=([ACTG]{25}GG[ACTG]{3}))', sequence)
matches_rev = re.finditer(r'(?=([ACTG]{25}GG[ACTG]{3}))', reverse_sequence)

# Load the appropriate model
if is_regression:
    rfm = dir + '/rfModelregressor.joblib'
    with open(rfm, 'rb') as f:
        rf = joblib.load(f)
else:
    rfm = dir + '/rfModelclassifier.joblib'
    with open(rfm, 'rb') as f:
        rf = joblib.load(f)
feature_lists = []
sequences = []


output_queue = Queue()
matches_queue = Queue(maxsize=num_threads*2)

print("Analysing")
with open(str(output_file), 'w') as f:
    f.write(header+str('\n'))

    for is_reverse, match_type in enumerate((matches, matches_rev)):
        is_reverse = bool(is_reverse)
        matches_process = Process(target=fill_queue, args=(match_type, matches_queue))
        matches_process.start()
        processes = [Process(target=score_sequences, args=(matches_queue, output_queue, is_reverse)) for x in range(num_threads)]
        for p in processes:
            p.start()
        num_done = 0
        while num_done < num_threads:
            output = output_queue.get()
            if output == "DONE":
                num_done += 1
            else:
                f.write(output)
                f.write("\n")
        matches_process.join()
        for process in processes:
            process.join()
