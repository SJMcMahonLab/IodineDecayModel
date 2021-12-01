# Calculate yield of SSBs based on a given set of physical and data files, and set energy and OH thresholds
from iodineCore import *

# Paths
folder = './' 							# Start from current folder, by default
codeOptions = ['opt2','opt4','opt6']	# One sub-folder for each physics model code

# Input file names
physFile = 'BaseDepositsOne.txt'
chemFile = 'BaseChemOne.txt'

# Physical threshold, chemistry interaction rate, and window size for SSB threshold
thresholds = [29.5, 31.75, 32.5]
chemProbs = [0.16, 0.16, 0.4]
delta = 0

# Data for scavenged treatment in Kandiya et al, Figure 3
bases=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]	
bases = list(map(int,bases))
fragmentProbs=[0.167971126,0.236717375,0.176224177,0.109325872,0.063830066,0.041236124,0.026909796,0.020033854,0.023767357,0.020595819,0.018215079,0.022275374,0.018536515,0.017420985,0.016371174,0.019224715]
logFragmentProbs = np.log(fragmentProbs)

def damageFileToFragments(depositArray,chemArray,bases,threshold, chemProb, delta=0):
	allBreaks = []
	lastBreaks = []
	lastlabel = -1
	anyBreak = 0
	for label, event in enumerate(zip(depositArray,chemArray)):
		eventDeposit, eventChem = event
		while label>lastlabel+1:
			allBreaks.append([])
			lastlabel+=1
		breaks = depositAndChemToBreaks(eventDeposit, eventChem, threshold, chemProb, delta)
		allBreaks.append(breaks)
		if len(breaks)>0:
			maxB = max(breaks)
			anyBreak+=1
		else:
			maxB = -1
		lastBreaks.append(maxB)

		lastlabel=label
	fragmentSizes, noBreaks = breaksToFragments(lastBreaks,len(bases))
	
	# If fragment sizes are too small, pad with minimum values
	while len(fragmentSizes)<len(bases): fragmentSizes = np.append(fragmentSizes,min(fragmentSizes))

	breakProb = anyBreak/len(lastBreaks)
	fullArray = [fragmentSizes[n] for n in bases]

	logArray = [np.log(fragmentSizes[n]) for n in bases]
	chiSqList = [pow(np.array(r)- np.array(f),2) for r,f in zip(logArray, logFragmentProbs)]
	chiSq = sum(chiSqList)

	return fullArray, chiSq

# Loop over models
for code, physThreshold, chemProb in zip(codeOptions,thresholds, chemProbs):
	fragments, chiSq = damageFileToFragments(depositArray,chemArray,bases,physThreshold, chemProb, delta)
	print(physThreshold,'\t',chemProb,'\t',chiSq,'\t\t','\t'.join(map(str,fullArray)) )
