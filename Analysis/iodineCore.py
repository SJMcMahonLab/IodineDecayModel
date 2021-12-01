# Core analysis tools for iodine energy deposit and break data. Handles phsp reading, manipulation, and some basic
# calcultions of relevant parameters and fragment size distributions. 

import os
import numpy as np
import pandas as pd
from collections import Counter

headers = ["MoleculeID or ParticlePDG", "Position X [um]", "Position Y [um]", "Position Z [um]", 
		   "EventID", "Track ID", "Step number", "Particle name", "Process name", "Volume name", 
		   "Volume copy number", "ParentA ID", "ParentB ID", "Vertex position x [um]", "Vertex position y [um]", 
		   "Vertex position z [um]", "Global time [ps]", "Energy deposited [keV]", "Kinetic energy [keV]"]

rng = np.random.default_rng(1)

verbose = False
cSize = 100000

#############################
# Basic data handling methods
#############################
def dataframeToBaseArray(baseArray,df):
	maxBase = df['Volume copy number'].max()
	maxEvent = df['EventID'].max()

	# Create first array
	if baseArray is None:
		baseArray = np.zeros( (maxEvent+1,maxBase+1) )
	else:
		# Pad array so it can hold all needed elements
		eventPadding = max(maxEvent+1-len(baseArray),   0)
		basePadding =  max(maxBase +1-len(baseArray[1]),0)
		baseArray = np.pad(baseArray,[(0,eventPadding) ,(0,basePadding)])

	for index, row in df.iterrows():
		baseArray[int(row['EventID'])][int(row['Volume copy number'])]+=1000*row['Energy deposited [keV]']

	return baseArray

def parseDataframe(filename, maxRows = None):
	baseArray = None
	df1 = pd.read_csv(filename, names = headers, delim_whitespace=True, nrows=maxRows, chunksize=cSize)
	for chunk in df1:
		baseArray = dataframeToBaseArray(baseArray,chunk)
	return baseArray

# Do same for chemistry. Base array size off physics array, so we need much less code here.
def parseDataframeToChem(filename, arraySize, maxTime = None, maxRows=None):
	baseArray = np.zeros(arraySize, dtype=np.int32)
	df1 = pd.read_csv(filename, names = headers, delim_whitespace=True, nrows=maxRows, chunksize=cSize)
	if maxTime == None: maxTime = np.inf
	for chunk in df1:
		for index, row in chunk.iterrows():
			# Filter out points with Z~0, to avoid radicals incident on back face of DNA
			if row['Position Z [um]']>1E-10 and row['Global time [ps]']<=maxTime:
				baseArray[int(row['EventID'])][int(row['Volume copy number'])]+=1
	return baseArray

def readChemistryDataframes(parentFolder, physFilename, chemFilename, maxTime = None, maxRows = None):
	baseArray = None
	chemArray = None
	listing = os.listdir(parentFolder)
	folders = sorted([f for f in listing if os.path.isdir(parentFolder+f)])
	print(folders)

	withFiles = []
	noFiles = []
	refShape = None
	for folder in folders:
		print(folder)
		chemSize = os.path.getsize(parentFolder+folder+'/'+chemFilename)
		if chemSize==0:	continue

		tempBaseArray = parseDataframe(parentFolder+folder+'/'+physFilename, maxRows=maxRows)
		if refShape is None:
			refShape = tempBaseArray.shape
		else:
			if refShape[0]<tempBaseArray.shape[0]: refShape[0]=tempBaseArray.shape[0]
			if refShape[1]<tempBaseArray.shape[1]: refShape[1]=tempBaseArray.shape[1]

		tempChemArray = parseDataframeToChem(parentFolder+folder+'/'+chemFilename, refShape, 
											 maxTime=maxTime, maxRows=maxRows)

		if baseArray is not None:
			baseArray = np.concatenate((baseArray, tempBaseArray))
		else:
			baseArray = tempBaseArray

		if chemArray is not None:
			chemArray = np.concatenate((chemArray, tempChemArray))
		else:
			chemArray = tempChemArray

	return baseArray, chemArray

def saveBaseArray(baseArray, name, fmt='%.18e'):
	np.savetxt(name,baseArray, fmt=fmt)
	return

def loadBaseArray(name):
	baseArray = np.loadtxt(name)
	return baseArray

def averageEnergyDeposits(depositArray, verbose=False):
	avgData = np.mean(depositArray, axis=0)
	if verbose:	print('\t'.join(map(str,avgData) ) )

	lowerCI = np.percentile(depositArray, 15.9, axis=0)
	median  = np.percentile(depositArray, 50, axis = 0)
	upperCI = np.percentile(depositArray, 84.1, axis=0)

	stdev = np.std(depositArray, axis=0)

	if verbose: print('\t','\t'.join(map(str,stdev)))
	return avgData, lowerCI, median, upperCI

#############################
# Simple break probabilities
#############################
# Convert energy deposits to breaks
def depositsToBreaks(eventData,thresh, delta = 0):
	breaks = []
	for base, deposit in enumerate(eventData):
		if delta==0 and deposit>thresh:
			breaks.append(base)
		else:
			excessE = deposit-(thresh-delta)
			if excessE>0:
				prob = excessE/(2*delta)
				if prob>1 or rng.random()<prob:
					breaks.append(base)

	return breaks

# Convert energy deposits and chemical events to breaks
def depositAndChemToBreaks(eventData, chemData, thresh, chemProb, delta = 0):
	breaks = []
	for base, deposit in enumerate(eventData):
		if delta==0 and deposit>thresh:
			breaks.append(base)
		else:
			excessE = deposit-(thresh-delta)
			if excessE>0:
				prob = excessE/(2*delta)
				if prob>1 or rng.random()<prob:
					breaks.append(base)

	for base, radicals in enumerate(chemData):
		if base not in breaks:
			if radicals >0:
				breakProb = 1-pow(1-chemProb,radicals)
				if rng.random()<breakProb:
					breaks.append(base)

	return breaks
	
# For a list of breaks, map them to probabilities
def breaksToProbs(breakList, strandSize):
	breakCounts = np.zeros(strandSize)
	for exposure in breakList:
		for base in exposure: 
			breakCounts[base]+=1

	#print(breakList)
	breakProbs = breakCounts/len(breakList)
	#print(breakProbs)

	return breakCounts, breakProbs

# Calculate probability that *any* base is broken following a decay
def anyBaseBreakProb(eventData, thresh):
	brokenBases = [1 if eventDeposit.max()>thresh else 0 for eventDeposit in eventData]
	return sum(brokenBases)/len(eventData)

# Convert break list to distribution of fragments (position of most distant break)
def breaksToFragments(breakList, maxBase=None, kandaiyaNormalise = False):
	fragmentCount = Counter(breakList)
	fragmentSizes = np.zeros(max(fragmentCount.keys())+1)
	for size, count in zip(fragmentCount.keys(), fragmentCount.values()):
		if size>=0: fragmentSizes[size]=count

	# Count half a break for empty bins to reflect uncertainty
	fragmentSizes = [max(f,0.5) for f in fragmentSizes]

	if maxBase is not None and maxBase<len(fragmentSizes):
		fragmentSizes = fragmentSizes[:maxBase]
	fragmentSizes = fragmentSizes/sum(fragmentSizes)

	if kandaiyaNormalise:
		reachProb = np.zeros(len(fragmentSizes))
		newFragmentSizes = np.zeros(len(fragmentSizes))
		newFragmentSizes[-1]=fragmentSizes[-1]
		reachProb[-1]=1
		for n in range(len(fragmentSizes)-2,-1,-1):
			reachProb[n]=reachProb[n+1]*(1-newFragmentSizes[n+1])
			newFragmentSizes[n]=fragmentSizes[n]/reachProb[n]

		fragmentSizes=newFragmentSizes

	if verbose:	print('\t'.join(map(str,fragmentSizes)))
	return fragmentSizes, fragmentCount[-1]/sum(fragmentCount.values())
