#!/usr/bin/python

#This program creates a feature matrix given a fasta (.fa) and corresponding bed (.bed) file
#Written by Daniel Reti, November 2016

#USAGE: ./this.program.py file.fa [file.bed]
#If no bed file, then score is not in present in matrix
#If both fa and bed provided, then only sequences with both score and appear in .fa will be in matrix
#TODO: May need to swap around score and sign in bed split - see TODO below

from __future__ import print_function

#determines gc content of given sequence
def gc(seq):
	l = len(seq)
	n = 0
	for r in seq:
		if r == 'C' or r == 'G':
			n+= 1
	return round(float(n)/l * 100, 2)  
		
#determines number of a given base in a given sequence
def content(seq, base):
	n = 0
	for r in seq:
		if r == base:
			n+=1
	return n

def di_content(seq, di):
	for i in di:
		count = 0
		for j in range(len(seq)-1):
			if (seq[j:j+2] == i):
				count+= 1
		print("{:4}".format(count), end='')
	print('  ', end='')

def pam(seq, di, bases):
	l = [0] * 16
	for i in di:
		multiply = bases.index(seq[24])
		offset = bases.index(seq[27])
		l[multiply*4 + offset] = 1
	print('     '.join(str(n) for n in l), end='  ')
		

#Nucleotide distribution across sequence
#Position 1 is 5' end, Position 21 is the N in the NGG PAM (3' end)
#Sequence structure: 5' 20merNGG 3'
#C7 means a C at position 7
def nucleotide(seq, bases):
	l = [0] * (len(seq)) * 4
	for i in range(len(seq)):
		count = 4 * i + bases.index(seq[i]) 
		l[count] = 1
	print('   '.join(str(n) for n in l), end='  ')
	
#Dinucleotide distribution accross sequence
#CT7 means C at 7 and T at 8 
#structure outlined above 
def dinucleotide(seq, di):
	#-1 is since a sequence of length N has N-1 dinucleotides
	l = [0] * ((len(seq))-1) * 16
	for i in range((len(seq))-1):
		l[16 * i + di.index(seq[i:i+2])] = 1
	print('    '.join(str(n) for n in l), end='')
		
import sys

#global variables

def main(args):
	d = {}
	bases = ['A', 'C', 'G', 'T']
	di = []
	# This can be adjusted if the sequence length changes down the track
	seqLength = 30
	for a in bases:
		for b in bases:
			di.append(str(a)+str(b))

	d = args[0]
			
	#reads bed file (arg2 in cmd line)
	if (len(args) > 1):
		g = open(args[1], 'r')
		for line in g:
			#TODO: may need to swap around sign and score as required
			#if score last
			#(ch, start, end, name, sign, score) = line.split()
			#if sign last
			(ch, start, end, name, score, sign) = line.split()
			if name in d:			
				d[name]['score'] = score		
		g.close()

	#view feature matrix
	#headers

	#with bed file (scores)
	if (len(args) > 1):
		LAYOUT = "{!s:25} {!s:6} {!s:2} {!s:2} {!s:2} {!s:2} {!s:12}"
		print(LAYOUT.format("Name", "GC_", "A", "C", "G", "T", "Activity"), end='')
	#without bed file (no scores)
	else:
		LAYOUT_1 = "{!s:30} {!s:6} {!s:2} {!s:2} {!s:2} {!s:2}"
		print(LAYOUT_1.format("Name", "GC_", "A", "C", "G", "T"), end='')

	# range is from 1 up to sequence length (+1 accounts for shift from zero)
	for i in range(1, (seqLength+1)):
		for j in bases:
			print("{:3}".format(str(j)+str(i)), end='  ')

	for i in range(1, (seqLength+1) -1):
		for j in di:
			print("{:4}".format(str(j)+str(i)), end='  ')

	for i in range(len(di)):
		print("{:4}".format(di[i]), end='')

	for i in range(len(di)):
		s = str(di[i][0]) + 'GG' + str(di[i][1])
		print("{:5}".format(s), end='  ')

	#new line after header
	print()

	#matrix data with bed file
	if (len(args) > 1):
		for a in d:
			if 'score' in d[a]:
				print(LAYOUT.format(
						a, 
						str(gc(d[a]['seq'])),
						content(d[a]['seq'], 'A'),
						content(d[a]['seq'], 'C'), 
						content(d[a]['seq'], 'G'),
						content(d[a]['seq'], 'T'),  
						d[a]['score']
				), end='  ')
				nucleotide(d[a]['seq'], bases)
				dinucleotide(d[a]['seq'], di)
				di_content(d[a]['seq'], di)
				pam(d[a]['seq'], di, bases)
				print()
	else:	
	#matrix data without bed file
		for a in d:
			print(LAYOUT_1.format(
					a, 
					str(gc(d[a]['seq'])),
					content(d[a]['seq'], 'A'),
					content(d[a]['seq'], 'C'), 
					content(d[a]['seq'], 'G'),
					content(d[a]['seq'], 'T'),  
			), end='  ')
			nucleotide(d[a]['seq'], bases)
			dinucleotide(d[a]['seq'], di)
			di_content(d[a]['seq'], di)
			pam(d[a]['seq'], di, bases)
			print()


if __name__ == "__main__":
	main(sys.argv)
