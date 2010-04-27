# -*- coding: utf-8 -*-
import sys
import nuto

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#Waehle Datenformat: 0 -int , 1 - short
format=0

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% start the real test file                                      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myStructure = nuto.StructureGrid(3) 

#Einlesen mit int
if format == 0:
	#create matrix
	myMatrix = nuto.IntFullMatrix(8,1,[155, 156, 162, 189, 209, 232, 241, 223])
	#myMatrix.Convert2short()

	#read matrix from file
	myInput = nuto.IntFullMatrix()
	myInput.ImportFromVtkASCIIFile("../../intNutoInputsmall")
	#myInput = nuto.ShortFullMatrix(1,8,[155, 156, 162, 189, 209, 232, 241, 223])

	myNumberOfRows = myInput.GetNumRows()
	myNumberOfColumns = myInput.GetNumColumns()

	myCompareInput = myInput.GetBlock(0,0,8,1)

	#calculate error
	if ((myNumberOfRows!=200000)):
		print '[' + myNumberOfRows+ '] : False Number of Rows.'
		error = True;

	if ((myNumberOfColumns!=1)):
		print '[' + myNumberOfColumns+ '] : False Number of Columns.'
		error = True;

	#if ((myMatrix-myCompareInput).Max()>1e-8) or ((myMatrix-myCompareInput).Max()<-1e-8):
	if (myMatrix-myCompareInput).Abs().Max()[0]>1e-8 :
		print '[' + system,sys.argv[0] + '] : myInputFile is not correct.'
		error = True;

	
#Einlesen mit short
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if format == 1:
	#create matrix
	myShortMatrix = nuto.IntFullMatrix(8,1,[155, 156, 162, 189, 209, 232, 241, 223])
	myShortMatrix.Convert2short()

	#read matrix from file
	myShortInput = nuto.ShortFullMatrix()
	myShortInput.ImportFromVtkASCIIFile("../../NutoInput")

	myShortNumberOfRows = myInput.GetNumRows()
	myShortNumberOfColumns = myInput.GetNumColumns()
	print  myShortInput.GetNumRows() 
	print myShortInput.GetValue(0,0)
	print myInput.GetValue(1,0)
	print myInput.GetValue(2,0)

	myShortCompareInput = myShortInput.GetBlock(0,0,8,1)

	#calculate error
	if ((myShortNumberOfRows!=840000)):
	    print '[' + myShortNumberOfRows+ '] : False Number of Rows.'
	    error = True;

	if ((myShortNumberOfColumns!=1)):
	    print '[' + myShortNumberOfColumns+ '] : False Number of Columns.'
	    error = True;

	#if ((myShortMatrix-myShortCompareInput).Max()>1e-8) or ((myShortMatrix-myShortCompareInput).Max()<-1e-8):
	#if ((myShortMatrix-myShortCompareInput).Abs().Max()[0]>1e-8) :
	#		 print '[' + system,sys.argv[0] + '] : myInputFile is not correct.'
	#		 error = True;

#Grid structure
myStructure = nuto.StructureGrid(3) 

myModul = nuto.DoubleFullMatrix(255,1) 
#myGridInput = nuto.StructureGrid.ImportFromVtkASCIIFile("../../intNutoInputsmall")
count = 0
for count in range(0,130):
	myModul.SetValue(count,1,0.)
for count in range(131,161):
	myModul.SetValue(count,1,8300.)
for count in range(162,254):
	myModul.SetValue(count,1,11500.)


