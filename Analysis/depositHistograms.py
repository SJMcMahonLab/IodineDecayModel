# Build a histogram of energy deposit event for a given phase space file

from iodineCore import *

folder = './' 							# Start from current folder, by default
codeOptions = ['opt2','opt4','opt6']	# One sub-folder for each physics model code

# Names of phase space files. Can chance 'One' to 'Two' for second strand
physFilename = 'IodineDepositStrandOne.phsp'

# Analaysis options. MaxRows sets number of phase space rows to analyze for efficiency. 'None' analyses all.
maxRows = None
doAnalysis = True

# For each code, iterate overall all sub-folders containing data, and aggregate into a single list for histogramming
for code in codeOptions:
	print('\n',code)
	depositArray = None
	chemArray = None
	listing = os.listdir(mainFolder+code)

	# Check each sub-folder for data
	folders = sorted([f for f in listing if os.path.isdir(mainFolder+code+'/'+f)])

	maindf = None
	for folder in folders:
		theFile = mainFolder+code+'/'+folder+'/'+physFile
		df = pd.read_csv(theFile, names = headers, delim_whitespace=True, nrows=maxRows)
		if maindf is None:
			maindf = df['Energy deposited [keV]']*1000
		else:
			maindf = pd.concat([maindf, df['Energy deposited [keV]']*1000])

	# Summary statistics
	if doAnalysis:
		print("Max E:\t",max(maindf))
		print("Total E:\t",sum(maindf))
		print("Events below 0.5 eV:\t",len(maindf[maindf<0.5]))
		print("Deposit below 0.5 eV:\t",sum(maindf[maindf<0.5]))
		print("Below 0.5 eV fraction:\t",sum(maindf[maindf<0.5])/sum(maindf))

	maindf = maindf[maindf>0.5]
	if doAnalysis:
		print('\nForPlot:')
		print("Max E:\t",max(maindf))
		print("Total E:\t",sum(maindf))
		print("Total Events:\t",len(maindf))

	# Actual histogram output
	hist,bins = np.histogram(maindf, bins=140, range=[0,35])
	print('\t'.join(map(str,bins)))
	print('\t'.join(map(str,hist)))
