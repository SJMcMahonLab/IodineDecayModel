# Calculate yields of DSBs for a given simulation

from iodineCore import *

folder = './' 							# Start from current folder, by default
codeOptions = ['opt2','opt4','opt6']	# One sub-folder for each physics model code
# Names of deposit array files, one for each strand
fileOne = 'BaseDepositsOne.txt'
fileTwo = 'BaseDepositsTwo.txt'

thresholds = [29.5, 31.75, 32.5] # Best fit energy thresholds

def depositDataToDSBs(strandOne, strandTwo, threshold, sep=10):
	totalDSBs = 0
	totalS1SSBs = 0
	totalS2SSBs = 0
	for event,(s1, s2) in enumerate(zip(strandOne, strandTwo)):
		s1Breaks = [n for n,d in enumerate(s1) if d>threshold ]
		s2Breaks = [n for n,d in enumerate(s2) if d>threshold ]
		ssbCount = len(s1Breaks)+len(s2Breaks)

		totalS1SSBs += len(s1Breaks)
		totalS2SSBs += len(s2Breaks)

		if event<5:
			print(s1Breaks, s2Breaks,ssbCount)

		remainingSSB = 0
		dsbList = []
		n=len(s1Breaks)-1
		while n>=0:
			ssb1 = s1Breaks[n]
			inDSB= False
			for m in range(len(s2Breaks)-1,-1,-1):
				ssb2=s2Breaks[m]
				if abs(ssb1-ssb2)<=sep:
					s2Breaks.pop(m)
					if inDSB:
						dsbList[-1][1].append(ssb2)
					else:
						inDSB=True
						dsbList.append([[ssb1],[ssb2]])

			s1Breaks.pop(n)
			n-=1
			if not inDSB:
				remainingSSB+=1

		remainingSSB+=len(s2Breaks)
		dsbCount = len(dsbList)

		totalDSBs+=len(dsbList)

		if event>=5:
			break

	return totalS1SSBs/strandOne.shape[0],totalS2SSBs/strandOne.shape[0], totalDSBs/strandOne.shape[0]


strandOneArrays = [loadBaseArray(folder+code+fileOne) for code in codeOptions]
strandTwoArrays = [loadBaseArray(folder+code+fileTwo) for code in codeOptions]
print('Opt2 SSB1\tOpt2 SSB2\tOpt2 DSB\tOpt4 SSB1\tOpt4 SSB2\tOpt4 DSB\tOpt6 SSB1\tOpt6 SSB2\tOpt6 DSB')
for n in range(3):
	s1SSB, s2SSB, dsbs = depositDataToDSBs(strandOneArrays[n], strandTwoArrays[n], threshold)
	print(s1SSB,'\t',s2SSB,'\t',dsbs,end='\t')


