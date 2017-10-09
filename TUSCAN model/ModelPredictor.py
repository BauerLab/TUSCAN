#!/usr/bin/python

#Given a sequence, this program uses a Random Forest to predict the activity 
#via a feature matrix, using a predefined master model created using a training set

# Can use a bed file, fasta file, or plain list of sequences

import sys
import FeatureMatrix
import pickle
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.externals import joblib
import numpy
import pybedtools
import argparse
from collections import OrderedDict

#Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', required=False, help='Genome')
parser.add_argument('-o', required=False, help='Output file')
parser.add_argument('-i', required=True, help='Input file')
parser.add_argument('-t', required=True, help='Type of input file (bed/txt/fa)')
parser.add_argument('-m', required=True, help='Type of Model: (Regression/Classification)')

args = parser.parse_args()

#assume file has name/ID included
has_name = True

# If a bed file is given
if (args.t == 'bed'):
	if not args.g:
		print('Genome required with bed file')
		sys.exit()
	b = open(args.i, 'r')
	bed_type = b.readline()
	bed_type = bed_type.split()
	bed_columns = len(bed_type)
	if bed_columns < 3:
		sys.stderr.write('Invalid bed file, must have at least 3 columns\n')
		sys.exit()
	c = pybedtools.BedTool(args.i)
	genome = args.g
	if bed_columns > 3:
		#Bed file contains Name column
		has_name = True
		d = c.sequence(name=True, fi=genome)
	else:
		#Bed file doesn't contain Name column
		has_name = False
		d = c.sequence(fi=genome)

#arguments to pass to matrix builder
name = str(args.i[:-len(str(args.t))-1]) + '_matrix.txt'

#get important features
l = []
if args.m == 'Regression':
	#f = open(args.f)
	f = open('rf_features_regression.txt')
elif args.m == 'Classification':
	f = open('rf_features_classification.txt')
else:
	sys.stderr.write('Invalid model type, must be Classification or Regression\n')
	sys.exit()	
for line in f:
	feat = line.split()
	for b in feat:
		l.append(b.strip('"'))
f.close()


#Valid DNA letters
valid = 'ACTG'
#Keeps input order the same
s = OrderedDict()
#Store sequences for recall later
if (args.t == 'bed' or args.t == 'fa'):
	if args.t == 'bed':
		f = open(d.seqfn, 'r')
	else:
		f = open(args.i, 'r')
	z = 1
	for line in f:
		line = line.rstrip()
		if line[0] == '>':
			if has_name:
				k = line[1:]
			else:
				k = (z+1)/2
		else:
			line = line.upper()
			if not all(i in valid for i in line):
				count = z
				if (args.t == 'bed'):
					count = z/2
				sys.stderr.write('The sequence at line ' + str(count) + ' contains an invalid nucleotide\n')
			elif len(line) != 30:
				count = z
				if (args.t == 'bed'):
					count = z/2
				sys.stderr.write('The sequence at line ' + str(count) + ' is not 30 base pairs, it is ' + str(len(line)) + '\n')
			elif line[25:27] != 'GG':
				revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
				reverseLine = revcompl(line)
				if reverseLine[25:27] != 'GG':
					count = z
					if (args.t == 'bed'):
						count = z/2
					sys.stderr.write('The sequence at line ' + str(count) + ' does not have a PAM motif\n')
				else:
					s[k] = {}
					s[k]['seq'] = reverseLine
					s[k]['dir'] = '-'
					count = z
					if (args.t == 'bed'):
						count = z/2
					sys.stderr.write('The sequence at line ' + str(count) + ' was in the negative orientation\n')
			else:
				s[k] = {}
				s[k]['seq'] = line
				s[k]['dir'] = '+'
		z += 1
	f.close()
elif (args.t == 'txt'):
	z = 1
	f = open(args.i, 'r')
	for line in f:
		line = line.rstrip()
		line = line.upper()
		if not all(i in valid for i in line):
			count = z
			if (args.t == 'bed'):
				count = z/2
			sys.stderr.write('The sequence at line ' + str(count) + ' contains an invalid nucleotide\n')
		elif len(line) != 30:
			sys.stderr.write('The sequence at line ' + str(z) + ' is not 30 base pairs, it is ' + str(len(line)) + '\n')
		elif line[25:27] != 'GG':
			revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])
			reverseLine = revcompl(line)
			if reverseLine[25:27] != 'GG':
				sys.stderr.write('The sequence at line ' + str(z) + ' does not have a PAM motif\n')
			else:
				s[z] = {}
				s[z]['seq'] = reverseLine
				s[z]['dir'] = '-'
				sys.stderr.write('The sequence at line ' + str(z) + ' was in the negative orientation\n')
		else:		
			s[z] = {}	
			s[z]['seq'] = line
			s[z]['dir'] = '+'
			z += 1
	f.close
else:
	print('Usage [-t bed/fa/txt]')
	sys.exit()

#Read in sequence data to create feature matrix
#generate feature matrix
stdout_ = sys.stdout
sys.stdout = open(name, 'w')
FeatureMatrix.main([s])
sys.stdout.close()
sys.stdout = stdout_

#grabs list of features
lines = []
with open(name, 'r') as f:
	features = f.readline()
features = features.split()
num_features = len(features)
features = features[1:]

#gets index of important features
a = [features.index(i) for i in l]

data = numpy.genfromtxt(name, dtype = 'f8', skip_header = 1, usecols = range(1, num_features))
train = data[:, a]

LAYOUT = '{!s:50} {!s:31} {!s:15} {!s:3}'
header = LAYOUT.format('ID', 'Sequence', 'Score', 'Dir')
#Open and predict on the randomForest
if args.m == 'Regression':
	with open('rfModelregressor.joblib', 'rb') as f:
		rf = joblib.load(f)
	scores = rf.predict(train)
	if args.o:
		with open(str(args.o), 'w') as f:
			f.write(header)
			for idx, a in enumerate(s):
				f.write(LAYOUT.format(a, s[a]['seq'], scores[idx], s[a]['dir']))
	else:
		print(header)
		for idx, a in enumerate(s):
			print(LAYOUT.format(a, s[a]['seq'], scores[idx], s[a]['dir']))


elif args.m == 'Classification':
	with open('rfModelclassifier.joblib', 'rb') as f:
		rf = joblib.load(f)
	scores = rf.predict(train)
	if args.o:
		with open(str(args.o), 'w') as f:
			f.write(header)
			for idx, a in enumerate(s):
				f.write(LAYOUT.format(a, s[a]['seq'], scores[idx], s[a]['dir']))
	else:
		print(header)
		for idx, a in enumerate(s):
			print(LAYOUT.format(a, s[a]['seq'], scores[idx], s[a]['dir']))


else:
	print('Usage: [-m Regression] or [-m Classification]')
