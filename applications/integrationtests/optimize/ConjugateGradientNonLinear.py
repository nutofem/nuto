import sys
import nuto
import os
import math
import numpy as np

# if set to true, the result will be generated (for later use in the test routine)
# otherwise, the current result will be compared to the stored result
createResult = False

# show the results on the screen
printResult = True

basePath = os.path.dirname(sys.argv[0])
# path in the original source directory and current filename at the and
pathToResultFiles = os.path.join(basePath, "results/")

# no error in file, modified, if error is detected
error = False


# start the real test file
def SetParameters(rParameterTuple):
    global x, y
    x = rParameterTuple[0]
    y = rParameterTuple[1]
    if (printResult):
        print ("[Python - Parameters routine]")
        print (x, y)


def Objective():
    o = scaleFactorX/12*(x-2)*(x-2)*(x-2)*(x-2)+0.5*(y-3)*(y-3)
    if (printResult):
       print ("[Python - Objective routine] (x,y)=(" + str(x) + "," + str(y) + ") : obj=" + str(o))
    return o


def Gradient():
    rGradient = [scaleFactorX/3*(x-2)*(x-2)*(x-2), (y-3)]
    if (printResult):
        print ("[Python - Gradient routine]")
        print (rGradient)
    return rGradient


def Hessian():
    rHessian = [scaleFactorX*(x-2)*(x-2), 0, 0, 1]
    if (printResult):
        print ("[Python - Hessian routine]")
        print (rHessian)
    return rHessian


x = 2.2
y = 15.0

scaleFactorX = 1e3

InitParameters = np.zeros((2, 1))
InitParameters[0, 0] = x
InitParameters[1, 0] = y

PythonCallback = nuto.CallbackHandlerPython()
PythonCallback.SetCallbackFunctions(SetParameters, Objective, Gradient, Hessian)

myOptimizer = nuto.ConjugateGradientNonLinear(2)
myOptimizer.SetCallback(PythonCallback)
myOptimizer.SetParameters(InitParameters)

myOptimizer.SetMaxFunctionCalls(1000000)
myOptimizer.SetMaxGradientCalls(1000)
myOptimizer.SetMaxHessianCalls(1000)
myOptimizer.SetMaxIterations(1000)
myOptimizer.SetMinObjective(0)
myOptimizer.SetMinDeltaObjBetweenRestarts(1e-6)
myOptimizer.SetAccuracyGradient(1e-12)

if (printResult):
    myOptimizer.SetVerboseLevel(10)
else:
    myOptimizer.SetVerboseLevel(0)

returnValue = myOptimizer.Optimize()
if (printResult):
    print ("Final objective : ", myOptimizer.GetObjective())
if (createResult):
    f = open(pathToResultFiles+'Objective.txt', 'w')
    f.write('#Correct Objective\n')
    f.write(str(myOptimizer.GetObjective()))
    f.close()
else:
    f = open(pathToResultFiles+'Objective.txt', 'r')
    f.readline()
    objectiveExact = float(f.readline())
    f.close()
    if (math.fabs(myOptimizer.GetObjective() - objectiveExact) > 1e-8):
        print ('[' + system,sys.argv[0] + '] : objective is not correct.')
        error = True

parameters = myOptimizer.GetParameters()
parameters = parameters.squeeze()
if (printResult):
    print ("Final Set of Parameters\n", parameters)
if createResult:
    np.savetxt(pathToResultFiles+"Parameters.txt", parameters, header="#Correct result")
else:
    ParametersExact = np.loadtxt(pathToResultFiles+"Parameters.txt", skiprows=1)
    print (ParametersExact)
    if (np.max(np.abs(ParametersExact - parameters)) > 1e-8):
        print ('[' + system, sys.argv[0] + '] : Parameters is not correct.')
        error = True

if (error):
    sys.exit(-1)
else:
    sys.exit(0)
