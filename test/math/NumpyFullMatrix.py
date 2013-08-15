import sys
import nuto
import ctypes
import numpy


myeigen= nuto.DoubleFullMatrix(2,2);

### take numpy convert to c++ and assign matrix to this numpy array and print matrix
A = numpy.array([[1,3.3] ,[0,1]],dtype='float64')
myeigen.convrtNumpyToMatrix(A)
#myeigen.printmatrix()

### access matrix in c++ here in python and print it
B =numpy.array([[1,1.0],[1,1.0]],dtype='float64')
myeigen.convrtMatrixToNumpy(B)
print B

## testing
c=A-B
d=numpy.array([[0,0] ,[0,0]],dtype='int32')
if numpy.array_equal(c, d):
 sys.exit(0)
else :
 sys.exit(-1)

