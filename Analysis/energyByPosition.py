# Plot energy deposit as a function of position along a strand

from iodineCore import *
import datetime

folder = './' 							# Start from current folder, by default
codeOptions = ['opt2','opt4','opt6']	# One sub-folder for each physics model code

# Names of deposit array files. Generate using dataFrameToArrays file on raw data phase spaces
filename = 'BaseDepositsOne.txt'

print('# Run at:\t',datetime.datetime.now())
print('# Folder:\t',folder)
print('# Options:\t',codeOptions)
for code in codeOptions:
	depositArray = loadBaseArray(folder+code+filename)
	print(code,end='\t')
	averageEnergyDeposits(depositArray, verbose=True)