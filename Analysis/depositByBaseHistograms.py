# Build a histogram of energy deposits in a single base

from iodineCore import *
import matplotlib.pyplot as plt
import datetime

folder = './' 							# Start from current folder, by default
codeOptions = ['opt2','opt4','opt6']	# One sub-folder for each physics model code

# Analysis file name, using by-base array
filename = 'BaseDepositsOne.txt'

# ID of base to analyse
baseID = 10

print('# Run at:\t',datetime.datetime.now())
print('# Folder:\t',folder)
print('# Options:\t',codeOptions)
for code in codeOptions:
	print(code)
	depositArray = loadBaseArray(folder+code+filename)
	baseData = depositArray[:,baseID]
	print('Events:\t',len(baseData))
	baseData = [b for b in baseData if b>0]
	print('Events with deposit:\t',len(baseData))
	print('Max Deposit:\t',max(baseData))
	print('Average Deposit:\t',np.mean(baseData))
	
	n,bins, patches = plt.hist(baseData, bins=400, range=(0,100))

	print('\t'.join(map(str,bins)))
	print('\t'.join(map(str,n)))
	print()
