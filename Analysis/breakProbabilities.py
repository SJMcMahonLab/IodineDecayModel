# Calculate break probabilities as a function of position for a given damage file

from iodineCore import *
import datetime

folder = './' 							# Start from current folder, by default
codeOptions = ['opt2','opt4','opt6']	# One sub-folder for each physics model code

# Names of deposit array files. Generate using dataFrameToArrays file on raw data phase spaces
filename = 'BaseDepositsOne.txt'

# Energy thresholds per break for each model respectively
thresholds = [29.5,31.75,32.50]

print('# Run at:\t',datetime.datetime.now())
print('# Folder:\t',folder)
print('# Options:\t',codeOptions)
print('# Thresholds:\t',thresholds)
for code,thresh in zip(codeOptions,thresholds):
	depositArray = loadBaseArray(folder+code+filename)
	print(code,end='\t')
	allBreaks = []
	for label, eventDeposit in enumerate(depositArray):
		breaks = depositsToBreaks(eventDeposit, thresh)
		allBreaks.append(breaks)

	rates = [len(b) for b in allBreaks]
	print(np.mean(rates),'\t',np.std(rates),'\t',end='\t')
	breakCounts, breakProbs = breaksToProbs(allBreaks, len(depositArray[0]))
	print('\t'.join(map(str,breakProbs)))
