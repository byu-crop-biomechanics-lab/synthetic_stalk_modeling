from part import *
from material import *
from section import *
from optimization import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

#define moduli values to be used in the for loop
Rind = 14600.0
Rindmoduli = [Rind, Rind*1.1, Rind*0.99]

#clear any old file with the same name and write header
fileloc = "\\auto\\fsi\\ryanal\\Crop Biomechanics\\Abaqus\\Results.txt"
FileResultsX=open(fileloc,'w')
FileResultsX.write("Modulus\t\tForce\n")
FileResultsX.close()

#loop for each material value
for rindmod in Rindmoduli:

	#assign material properties
	mdb.models['Model-1'].materials['Rind'].elastic.setValues(table=((rindmod, 0.3), ))


	#create job
	myJob = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
		explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
		memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
		multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
		numCpus=1, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
		userSubroutine='', waitHours=0, waitMinutes=0)
	
	
	#run job and wait to completion
	mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
	myJob.waitForCompletion()

	#query reaction force from output database
	odb = openOdb(path='Job-1.odb')
	numFrame=odb.steps['Step-2'].frames[-1]
	RForce=numFrame.fieldOutputs['RF']
	regS1 = odb.rootAssembly.nodeSets['READSET']

	FX = RForce.getSubset(region=regS1).values[0].data[2]
	
	#write to file (in append mode)
	FileResultsX=open(fileloc,'a')
	FileResultsX.write('%.2E\t' % (rindmod))
	FileResultsX.write('%.2E\t' % (FX))
	FileResultsX.write('\n')
	
	#close the file and output database, and then delete the output database
	FileResultsX.close()
	odb.close()
	del mdb.jobs['Job-1']
