import sys
import nuto
import os
import math
import ctypes
import numpy


myeigen= nuto.IntFullMatrix(2,2);

### take numpy convert to c++ and assign matrix to this numpy array and print matrix
A = numpy.array([[1,2] ,[0,1]],dtype='int32')
myeigen.convrtfn(A)
myeigen.printmatrix()

### access matrix in c++ here in python and print it
B =numpy.array([[1,1],[1,1]],dtype='int32')
#myeigen.convrtToPy(B)
print B

## testing
c=A-B
#d=numpy.array([[0,0] ,[0,0]],dtype='int32')
#if numpy.array_equal(c, d):
# sys.exit(0)
#else :
 #sys.exit(-1)

