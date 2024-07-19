# ############################################################################
#
# Code to parse a simple event, position, DSB file into an SDD format.
#
# ############################################################################
import sys
import os

import chromModel
chromosomes = 46
radius = 4.65 # Microns
chromModel.subDivideSphere(chromosomes, radius)

# Write SDD header
def writeHeader(f, energy, incident, dose, geometry, DNADensity, bdRange, damageCount, primaryCount):
	sourceType = 1
	f.write("SDD version, SDDv1.0;\n")
	f.write("Software, Model Name;\n")
	f.write("Author, Name, Email, Data;\n")
	f.write("Simulation Details, ;\n")
	f.write("Source, ;\n")
	f.write("Source type, "+str(sourceType)+";\n")
	if type(incident) is list: incident = ','.join(map(str,incident))
	f.write("Incident particles, " +str(incident)+";\n")
	if type(energy) is list: energy = ','.join(map(str,energy))
	f.write("Mean particle energy, "+str(energy)+";\n")
	f.write("Energy distribution, M, "+str(0)+";\n")
	f.write("Particle fraction, "+str(1.0)+";\n")
	f.write("Dose or fluence, "+str(dose)+", 1;\n")
	f.write("Dose rate, 0.0;\n")
	f.write("Irradiation target, Simple spherical cell model, radius "+str(radius)+" um;\n")
	f.write("Volumes, 0,5,5,5,0,0,0, 1,"+str(radius)+","+str(radius)+","+str(radius)+","+str(radius)+",0,0,0;\n")
	f.write("Chromosome sizes, 46, "+ ", ".join(map(str,[6.1E3/46 for n in range(46)])) +";\n")
	f.write("DNA Density, "+str(DNADensity)+";\n")
	f.write("Cell Cycle Phase, 0;\n")
	f.write("DNA Structure, 0, 1;\n")
	f.write("In vitro / in vivo, 0;\n")
	f.write("Proliferation status, 1,;\n")
	f.write("Microenvironment, 20, 0.01;\n")
	f.write("Damage definition, 0, 0, 10, "+str(bdRange)+", 0, 0;\n")
	f.write("Time, 0;\n")
	f.write("Damage and primary count, "+str(damageCount)+", "+str(primaryCount)+";\n")
	f.write("Data entries, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;\n")
	f.write("Additional information,;\n")
	f.write("***EndOfHeader***;\n\n")

# Write hits to file
def writeHits(f, hits):
	for hitSet in hits:
		for event in hitSet:
			f.write('; '.join(map(str,event)))
			f.write(';\n')

# Write dataset to file
def writeToFile(hits, outFile, target='Unspecified', geometry=[0,1,1,1],
				DNADensity=14.484, bdRange=-1, incident=0, energy=2400.0, dose=1.0):
	damageCount = sum([len(h) for h in hits])
	primaryCount = sum([sum([min(1,int(event[0][0])) for event in hitSet]) for hitSet in hits])
	with open(outFile,'w') as f:
	 	writeHeader(f, energy, incident, dose, geometry, DNADensity,
	 				bdRange, damageCount, primaryCount)
	 	writeHits(f, hits)

# ##############################
# Pick folder to work on here
# ##############################
#folder = 'FolderPath'

def convertInFolder():
	files = list(os.listdir())
	files = [f for f in files if (f[-4:]=='.txt' and 'sdd' not in f)]

	for inFile in files:
		outName = inFile.rsplit('.',1)[0]+'_sdd.txt'
		data = []
		eventNo = 0
		with open(inFile) as f:
			chromModel.subDivideSphere(chromosomes, radius) # Randomise chromosomes by file
			for row in f:
				if int(row[0])>0: data.append([])
				rowData = row.strip().split()

				particle = str(int(rowData[0]))+','+str(eventNo)
				position = ','.join(rowData[1:4])
				chromSpec = chromModel.modelChromosome(*list(map(float,rowData[1:4])))
				spec = '0, '+str(int(rowData[4]))+', 1'
				data[-1].append([particle, position, chromSpec[0], chromSpec[1], spec])

		writeToFile(data,outName)


convertInFolder() # Can add folder name argument to work in a different folder