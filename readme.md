
# Iodine-125 DNA damage analysis tools

This is a set of programs used in the simulation and analysis of energy deposition and DNA strand breaking, based on published data by (Kandaiya 1996) and modelling approches by (Nikjoo 1996). These models are re-implemented in Topas-nBio, and used to simulated energy depositions with a range of physical models of electron energy deposition.

This code has two components - a Topas extension and parameter files used to generate physical damage simulations, and a set of python code used to analyse the output of these simulations. These parts are briefly described below.

## Topas code and parameter files

To enable scoring of radical interactions in the same manner as (Nikjoo 1996), a scoring extension, NikjooTuple, is provided. This scores the first interaction of OH radicals with backbone strands, and kills all radicals interacting with DNA to simulate scavenging. This must be compiled with Topas as an extension to use the model parameter files.

A set of parameter files are also provided. **FullIodineChemistry.txt** is the main file used to simulate interactions in this work. Running it simulates 125I decays in a DNA strand, logging energy depositions and chemical interactions throughout the strand, and scoring these to phase spaces for subsequent anlysis. Include files chemistryTuple.txt and physicsTuple.txt define their respective scorers, enabling them to be readily turned off if desired. **TOPASChemistry.txt** defines standard rate constants for relevant chemical processes.

The example file is set to simulate with Geant4-DNA physics opt2, but this can be changed to opt4 and opt6 to explore the impact of other models.

It should be noted that in some versions of Topas there can be rare crashes when tracking chemical species with scavenging effects included. These are uncommon and not thought to significantly impact on simulation results, but may make it advisable to run a large number of shorter simulations to provide good statistics.

## Analysis code

Analysis code is provided, written in python3. This takes the output of the Topas model and generates relevant biological parameters, such as the yield of strand breaking. It depends on the following packages:

- numpy
- pandas
- matplotlib

The main analysis methods are stored in the iodineCore.py file, which includes methods to process phase space files, generate distributions of strand breaks and DNA fragments, and format some of this data output. In addition to be usable in new code, a number of reference files are provided to replicate key analyses in paper. These are:

- **dataFrameToArrays.py**: This is a basic data handling package, which converts folders of raw phase space data into per-event, per-base arrays to facilitate other subsequent analyses and significantly improve performance. Example file paths are provided, for simulations run with Geant4-DNA models opt4 and opt6, on the assumption that individual runs are grouped into subfolders of the form "`./opt2/folder1/`", "`./opt2/folder2/`" etc. and similarly for other physics lists. Data is collected from all subfolders and combined into a single output file for each physics process. 

- **energyByPosition.py**: This is a basic data analysis tool, which can be used to calculate the average energy deposit (or number of free radical interactions) per base in a given processed data array.

- **ssbChem.py**: This file calculates the yield of SSB for a given set of physical and chemical interaction files, based on provided break energy thresholds and OH interaction rates. This is generated as a distribution of fragment sizes for each model, together with the sum-square difference with the Kandaiya et al reference data.

- **breakProbabilitics.py**: This calculates the raw probability of a physical SSB in a DNA strand based on a single energy threshold.

- **depositHistograms.py** and **depositByBaseHistograms.py**: These tools generate histograms of energy deposits, either on a per-interaction basis (depositHistograms) or the amount of energy deposited in a given base on a per-event basis (depositByBaseHistograms.py).

- **dsbModel.py**: This tool calculates rates of DSBs, using as input energy depositions within both strands of the DNA. DSBs are called if two base pairs on opposite strands differ by a distance of fewer than 'sep' base pairs (default 10, but can be set as an argument).

## Contact

For questions/comments/bug reports, please contact stephen.mcmahon (at) qub.ac.uk

## References

[Kandaiya 1996] Kandaiya S, Lobachevsky PN, D’Cunha G, Martin RF. DNA strand breakage by 125I-decay in a synthetic oligodeoxynucleotide. Fragment distribution and evaluation of DMSO protection effect. Acta Oncol (Madr). 1996;35:803–8. 

[Nikjoo 1996] Nikjoo H, Martin RF, Charlton DE, Terrissol M, Kandaiya S, Lobachevsky P. Modelling of Auger-induced DNA damage by incorporated 125I. Acta Oncol (Madr). 1996;35:849–56. 