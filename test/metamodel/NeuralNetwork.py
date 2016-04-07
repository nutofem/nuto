import nuto
import random
import math
import tempfile
import sys
import os

#if set to true, the result will be generated (for later use in the test routine)
#otherwise, the current result will be compared to the stored result
#createResult = True
createResult = False

#show the results on the screen
printResult = True

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
#%  Random seed for training data generation and initial weights %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randomSeed = 1234567
random.seed(randomSeed)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%  Generate Training data in the range min max                  %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minCoordinate = -3
maxCoordinate =  3
meanNoise = 0
sigmaNoise = 1e-1
dimInput = 2
dimOutput = 2
numSamples = 10
numNeuronsHiddenLayer = 3
verboseLevel = 0

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% some auxiliary routines                                       %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def RandomSamplesUniform(minCoordinate,maxCoordinate,dim,numSamples):
    sample = nuto.DoubleFullMatrix(dim,numSamples)
    for count1 in range(0, dim):
        for count2 in range(0, numSamples):
            sample.SetValue(count1,count2,random.uniform(minCoordinate,maxCoordinate))
    return sample

def RandomSamplesGauss(mean,sigma,dim,numSamples):
    sample = nuto.DoubleFullMatrix(dim,numSamples)
    for count1 in range(0, dim):
        for count2 in range(0, numSamples):
            sample.SetValue(count1,count2,random.gauss(mean,sigma))
    return sample

def ExactFunction(dimOutput, coordinates):
    numSamples = coordinates.GetNumColumns()
    sample = nuto.DoubleFullMatrix(dimOutput,numSamples)
    for countdim in range(0, dimOutput):
        for count in range(0, numSamples):
            #sample.SetValue(countdim,count,(coordinates.GetValue(countdim,count)*coordinates.GetValue((countdim+1)%dimInput,count)))
            sample.SetValue(countdim,count,(coordinates.GetValue((countdim+1)%dimInput,count)))
    return sample

def wait(str=None, prompt='Press return to continue...\n'):
    raw_input(prompt)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% create support points                                         %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SupportPointsInput = RandomSamplesUniform(minCoordinate,maxCoordinate,dimInput,numSamples)
SupportPointsOutputExact = ExactFunction(dimOutput,SupportPointsInput)
SupportPointsOutputNoise = RandomSamplesGauss(meanNoise,sigmaNoise,dimOutput,numSamples)
SupportPointsOutputWithNoise = SupportPointsOutputExact + SupportPointsOutputNoise

#create network with certain number of neurons in each hiddenlayer
hiddenLayerNumNeurons = nuto.IntFullMatrix(1,1)
hiddenLayerNumNeurons.SetValue(0,0,numNeuronsHiddenLayer)
myNetwork = nuto.NeuralNetwork(hiddenLayerNumNeurons)
myNetwork.InitRandomNumberGenerator(randomSeed)
#myNetwork.UnsetBayesianTraining()
myNetwork.SetBayesianTraining()
myNetwork.SetInitAlpha(1e-5)
myNetwork.SetAccuracyGradient(1e-5)
myNetwork.SetMinDeltaObjectiveBetweenRestarts(1e-4)
myNetwork.SetMinDeltaObjectiveBayesianIteration(1e-3)
myNetwork.SetShowSteps(1)
myNetwork.SetMaxFunctionCalls(0)  #this is not very practical, but only used for the test file
myNetwork.SetVerboseLevel(verboseLevel)
myNetwork.SetMaxBayesianIterations(1)
#myNetwork.UseFullHessian()  #forpreconditioning in the conjugate gradient method
myNetwork.UseDiagHessian()

myNetwork.SetSupportPoints(dimInput,dimOutput,SupportPointsInput,SupportPointsOutputWithNoise)

#set transformations for input/output
lowerBound = -1
upperBound = 1
myNetwork.AppendMinMaxTransformationInput(lowerBound,upperBound)
myNetwork.AppendMinMaxTransformationOutput(lowerBound,upperBound)
myNetwork.BuildTransformation()
#set transfer functions
myNetwork.SetTransferFunction(0,nuto.NeuralNetwork.TanSig)
myNetwork.SetTransferFunction(1,nuto.NeuralNetwork.PureLin)

#build/train the network
#pdb.set_trace()
myNetwork.Build()

if (printResult):
    Parameters = nuto.DoubleFullMatrix(1,1)
    myNetwork.GetParameters(Parameters)
    print 'final set of parameters:'
    Parameters.Trans().Info()

#calculate error for the training samples including standard deviation of approximation
SupportPointsApproximation = nuto.DoubleFullMatrix(1,1)
SupportPointsApproximationMin = nuto.DoubleFullMatrix(1,1)
SupportPointsApproximationMax = nuto.DoubleFullMatrix(1,1)
myNetwork.SolveConfidenceInterval(SupportPointsInput,SupportPointsApproximation,SupportPointsApproximationMin,SupportPointsApproximationMax)

#objective
objective = myNetwork.Objective()
if (printResult):
    print 'objective including regularization terms (transformed space):\n  ' + str(objective)
if (createResult):
    f = open(pathToResultFiles+'Objective.txt','w')
    f.write('#Correct Objective\n')
    f.write(str(objective))
    f.close()
else:
    f = open(pathToResultFiles+'Objective.txt','r')
    f.readline()
    objectiveExact = float(f.readline())
    f.close()
    if (math.fabs(objective-objectiveExact)>1e-8):
        print '[' + system,sys.argv[0] + '] : objective is not correct.(' + str(math.fabs(objective-objectiveExact)) + ')'
        print 'objectiveExact including regularization terms (transformed space):\n  ' + str(objectiveExact)
        error = True;
	
	
#gradient
NumParameters = myNetwork.GetNumParameters()
gradient = nuto.DoubleFullMatrix(NumParameters,1)
myNetwork.Gradient(gradient)
if (printResult):
    print 'gradient:'
    gradient.Trans().Info()
if (createResult):
    gradient.WriteToFile(pathToResultFiles+"Gradient.txt"," ","#Correct gradient matrix","  ")
else:
    gradientExact = nuto.DoubleFullMatrix(NumParameters,1)
    gradientExact.ReadFromFile(pathToResultFiles+"Gradient.txt",1," ")
    if ((gradientExact-gradient).Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : gradient is not correct.(' + str((gradientExact-gradient).Abs().Max()) + ')'
        print 'gradientExact:'
        gradientExact.Trans().Info()
        error = True;

#hessian
hessian = nuto.DoubleFullMatrix(NumParameters,NumParameters)
myNetwork.HessianFull(hessian)
if (printResult):
    print 'hessian:'
    hessian.Info()
if (createResult):
    hessian.WriteToFile(pathToResultFiles+"Hessian.txt"," ","#Correct hessian matrix","  ")
else:
    hessianExact = nuto.DoubleFullMatrix(NumParameters,1)
    hessianExact.ReadFromFile(pathToResultFiles+"Hessian.txt",1," ")
    RelError = (hessianExact-hessian).ElementwiseDiv(hessianExact)
    if (RelError.Abs().Max()>1e-8):
        print '[' + system,sys.argv[0] + '] : hessian is not correct.(' + str(RelError.Abs().Max()) + ')'
        print 'hessianExact:'
        hessianExact.Info()
        print 'relative Error:'
        RelError.Info(12,3)
        error = True;




if (error):
    print "This test fails because the random number generator in NuTo::Metamodel was changed from a mersenne twister implementation from an external lib to the c++11 version in lib<random>. Maybe with different parameters/seed/idk. However, the results are in the same order of magnitude. \nI suggest updating the reference files. \nOver and Out. TT\n"
    sys.exit(-1)
else:
    sys.exit(0)
