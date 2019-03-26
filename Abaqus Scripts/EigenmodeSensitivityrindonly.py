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


#set baseline rind modulus
stalknumber = '1_2min95'
E_rind = 14673.0



#set model names: first run baseline model (BASE), then vary elastic moduli (E1 & E3) and shear moduli (G12 & G13) for the rind (R)
models = ['BASE', 'R_E1', 'R_E3', 'R_G12', 'R_G13']


#clear any old file with the same name and write header
filename = stalknumber + '_Results.txt'
FileResultsX=open(filename,'w')
FileResultsX.write("Model\t\tLoad\n")
FileResultsX.close()


#loop through the different models for the sensitivity study
for idx, model in enumerate(models):
	
	#print the current model to the Abaqus script window
	print(' ')
	print('Now running model...')
	print(model)
	
	#initialize the rind and pith modulus values
	rind = [E_rind/100, E_rind, E_rind/1000, E_rind/10]
	# pith = [x / 100 for x in rind]
	
	#depending on the iteration, change the appropriate parameter by 1%
	if idx > 0 and idx < 5:
		rind[idx - 1] = rind[idx - 1]*1.01
		
	# if idx > 4:
	# 	pith[idx - 5] = pith[idx - 5]*1.01
	
	#print the modulus values to the Abaqus script window
	print('with material properties of:')
	print(rind)
	# print(pith)
	
	#assign the material properties
	mdb.models['Model-1'].materials['Rind'].elastic.setValues(table=((rind[0], rind[0], 
		rind[1], 0.0, 0.0, 0.0, rind[2], rind[3], rind[3]), ), type=
		ENGINEERING_CONSTANTS)
	
	# mdb.models['Model-1'].materials['Pith'].elastic.setValues(table=((pith[0], pith[0], 
	# 	pith[1], 0.0, 0.0, 0.0, pith[2], pith[3], pith[3]), ), type=
	# 	ENGINEERING_CONSTANTS)
	
	#set job name and output dat file name
	jobname = stalknumber + '-' + model
	datname = jobname + '.dat'
	
	#create job
	myJob = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
		explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
		memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
		multiprocessingMode=DEFAULT, name=jobname, nodalOutputPrecision=SINGLE, 
		numCpus=1, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
		userSubroutine='', waitHours=0, waitMinutes=0)
	
	#run job and wait to completion
	mdb.jobs[jobname].submit(consistencyChecking=OFF)
	myJob.waitForCompletion()

	#go through resulting *.dat text file and find eigenvalue readout
	myOutdf = open(datname,'r')
	stline=' MODE NO      EIGENVALUE\n'
	lines = myOutdf.readlines()
	ss=0

	for i in range(len(lines)-1):
		if lines[i] == stline :
			#print lines[i]
			ss=i

	f1=lines[ss+3]
	MinEigen=float(f1[15:24])
	myOutdf.close()

	#print eigenvalue to Abaqus script window
	print(MinEigen)
	
	#write to file (in append mode), and then close the file
	FileResultsX=open(filename,'a')
	FileResultsX.write(model)
	FileResultsX.write('\t\t%.5E\n' % (MinEigen))
	FileResultsX.close()
