import sys, string, os
from ROOT import TFile, TNtuple
from ROOT import TCanvas, TGraphErrors, TLegend
from ROOT import gROOT
from array import array
ifn = os.path.expandvars("FinalShallow.csv")
ofn = 'FinalShallow.root'
print ('opening file'), ifn, '...'
infile = open( ifn, 'r' )
lines  = infile.read().splitlines()
title  = "Mean Detected Photons Per Muon"
print ('writing file'), ofn, '...'
outfile = TFile( ofn, 'RECREATE', 'ROOT file with an NTuple' )
i=0
Arrx = {}
Det = {}
Errs = {}
Err = {}
labels = []
for line in lines:
	words = line.split(",")
	if len(words) == 1:
		labels.append(words)
		i+=1
		if i > 1:
			print("Would you like to make a histogram of %s? [y/n]" % (labels[i-2]))
			yes = ('yes','y','ye')
			no = ('no','n')
			name = labels[i-2]
			hist = str(labels[i-2])
			choice = raw_input()
			if choice in yes:
				name = TGraphErrors(len(Arrx[i-1]), Arrx[i-1], Det[i-1], Err[i-1], Errs[i-1])
				name.Write("%s" % (hist))
			elif choice in no:
				p = True
			else:
				sys.stdout.write("Please respond with 'yes' or 'no'")
		Arrx[i] = array( 'f', [])
		Det[i] = array( 'f', [])
		Errs[i] = array( 'f', [])
		Err[i] = array( 'f', [])
	else:
		Arrx[i].append(float(words[0]))
		Det[i].append(float(words[2]))
		Errs[i].append(float(words[3]))
		Err[i].append(float(words[1]))
print ('done')