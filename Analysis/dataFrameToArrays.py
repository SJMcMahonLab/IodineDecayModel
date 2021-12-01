# Convert raw phase space data frames to energy deposit and OH interaction arrays

from iodineCore import *

folder = './' 							# Start from current folder, by default
codeOptions = ['opt2','opt4','opt6']	# One sub-folder for each physics model code

# Names of phase space files. Can chance 'One' to 'Two' for second strand
physFilename = 'IodineDepositStrandOne.phsp'
chemFilename = 'IodineDepositChemOne.phsp'

baseOutname = 'BaseDepositsOne.txt'
chemOutname = 'BaseChemOne.txt'

# For each model code, aggregate phase space data into per-event arrays
for code in codeOptions:
	print(folder+code)
	baseArray, chemArray = readChemistryDataframes(folder+code+'/', physFilename, chemFilename, maxTime=1000.0)
	saveBaseArray(baseArray, code+baseOutname)
	saveBaseArray(chemArray, code+chemOutname, fmt="%2d")
