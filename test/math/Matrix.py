# -*- coding: utf-8 -*-
import sys
import nuto
import os
import math

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
#createResult = True
createResult = False

#show the results on the screen
printResult = False

#system name and processor
system = sys.argv[1]+sys.argv[2]

#path in the original source directory and current filename at the and
pathToResultFiles = os.path.join(sys.argv[3],"results",system,os.path.basename(sys.argv[0]))

#remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt,'')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% no error in file, modified, if error is detected              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = False

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% start the real test file                                      %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = nuto.DoubleVector(6)
v=[1.2,2,3,4,5,6]
a=nuto.DoubleFullMatrix(2,3,v)
if (printResult):
   print "a"
   a.Info(5)  #width, precision
   print ""
v2 = nuto.DoubleVector(6)
v2=[11.,12.,13.,14.,15.,16.]
b=nuto.DoubleFullMatrix(3,2,v2)
if (printResult):
    print "b"
    print b
    b.Info(5)  #width, precision
    print ""

result1 = a*b
if (printResult):
    print "result1 = a*b"
    result1.Info(5)
    print ""
if createResult:
    print pathToResultFiles
    result1.WriteToFile(pathToResultFiles+"Result1.txt"," ","#Correct result","  ")
else:
    result1Exact = nuto.DoubleFullMatrix(1,1)
    print pathToResultFiles
    result1Exact.ReadFromFile(pathToResultFiles+"Result1.txt",1," ")
    if ((result1Exact-result1).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result1 is not correct.'
        error = True;

result2 = a*a.Trans()
if (printResult):
    print "result2 = a*at"
    result2.Info(5)
    print ""
if createResult:
	result2.WriteToFile(pathToResultFiles+"Result2.txt"," ","#Correct result","  ")
else:
    result2Exact = nuto.DoubleFullMatrix(1,1)
    result2Exact.ReadFromFile(pathToResultFiles+"Result2.txt",1," ")
    if ((result2Exact-result2).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result2 is not correct.'
        error = True;

result3 = a.Trans()*a
if (printResult):
    print "result3 = at*a"
    result3.Info(5)
    print ""
if createResult:
	result3.WriteToFile(pathToResultFiles+"Result3.txt"," ","#Correct result","  ")
else:
    result3Exact = nuto.DoubleFullMatrix(1,1)
    result3Exact.ReadFromFile(pathToResultFiles+"Result3.txt",1," ")
    if ((result3Exact-result3).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result3 is not correct.'
        error = True;

result4 = a.Trans()*b.Trans()
if (printResult):
    print "result4 = at*bt"
    result4.Info(5)
    print ""
if createResult:
	result4.WriteToFile(pathToResultFiles+"Result4.txt"," ","#Correct result","  ")
else:
    result4Exact = nuto.DoubleFullMatrix(1,1)
    result4Exact.ReadFromFile(pathToResultFiles+"Result4.txt",1," ")
    if ((result4Exact-result4).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result4 is not correct.'
        error = True;

result5 = result3+result4
if (printResult):
    print "result5 = result3+result4"
    result5.Info(5)
    print ""
if createResult:
	result5.WriteToFile(pathToResultFiles+"Result5.txt"," ","#Correct result","  ")
else:
    result5Exact = nuto.DoubleFullMatrix(1,1)
    result5Exact.ReadFromFile(pathToResultFiles+"Result5.txt",1," ")
    if ((result5Exact-result5).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result5 is not correct.'
        error = True;

int_result6=result3.Convert2int()
if (printResult):
    print "int_result6 = convert result3 to integer "
    int_result6.Info(5)
    print ""
if createResult:
	int_result6.WriteToFile(pathToResultFiles+"Result6.txt"," ","#Correct result","  ")
else:
    int_result6Exact = nuto.IntFullMatrix(1,1)
    int_result6Exact.ReadFromFile(pathToResultFiles+"Result6.txt",1," ")
    if ((int_result6Exact-int_result6).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result6 is not correct.'
        error = True;

int_result7=result4.Convert2int()
if (printResult):
    print "int_result7 = convert result4 to integer "
    int_result7.Info(5)
    print ""
if createResult:
	int_result7.WriteToFile(pathToResultFiles+"Result7.txt"," ","#Correct result","  ")
else:
    int_result7Exact = nuto.IntFullMatrix(1,1)
    int_result7Exact.ReadFromFile(pathToResultFiles+"Result7.txt",1," ")
    if ((int_result7Exact-int_result7).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result7 is not correct.'
        error = True;

int_result8 = int_result6 + int_result7
if (printResult):
    print "int_result8 = int_result6 + int_result7"
    int_result8.Info(5)
    print ""
if createResult:
	int_result8.WriteToFile(pathToResultFiles+"Result8.txt"," ","#Correct result","  ")
else:
    int_result8Exact = nuto.IntFullMatrix(1,1)
    int_result8Exact.ReadFromFile(pathToResultFiles+"Result8.txt",1," ")
    if ((int_result8Exact-int_result8).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result8 is not correct.'
        error = True;

result9=int_result8.Convert2double()
if (printResult):
    print "result9 = convert int_result8 back to double"
    result9.Info(5)
    print ""
if createResult:
	result9.WriteToFile(pathToResultFiles+"Result9.txt"," ","#Correct result","  ")
else:
    result9Exact = nuto.DoubleFullMatrix(1,1)
    result9Exact.ReadFromFile(pathToResultFiles+"Result9.txt",1," ")
    if ((result9Exact-result9).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result9 is not correct.'
        error = True;

result10=result9*0.5
if (printResult):
    print "result10 = result9 scaled with factor 0.5"
    result10.Info(5)
    print ""
if createResult:
	result10.WriteToFile(pathToResultFiles+"Result10.txt"," ","#Correct result","  ")
else:
    result10Exact = nuto.DoubleFullMatrix(1,1)
    result10Exact.ReadFromFile(pathToResultFiles+"Result10.txt",1," ")
    if ((result10Exact-result10).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : result10 is not correct.'
        error = True;

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
