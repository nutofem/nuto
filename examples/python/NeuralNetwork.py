import nuto
import random
import math
#import numpy
#import Gnuplot
import tempfile
import sys
import os

#show the results on the screen
printResult = True

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% define correct results                                        %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#objectiveExact = 0
#gradientExact = nuto.DoubleFullMatrix()
#hessianExact = nuto.DoubleFullMatrix()
#gradientExact.restore(path_test_sources+"NeuralNetworkGradient_"+system+".txt",nuto.XML)

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
myNetwork.SetShowSteps(100)
myNetwork.SetMaxFunctionCalls(1000000)
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

#calculate error for the training samples including standard deviation of approximation
SupportPointsApproximation = nuto.DoubleFullMatrix(1,1)
SupportPointsApproximationMin = nuto.DoubleFullMatrix(1,1)
SupportPointsApproximationMax = nuto.DoubleFullMatrix(1,1)
myNetwork.SolveConfidenceInterval(SupportPointsInput,SupportPointsApproximation,SupportPointsApproximationMin,SupportPointsApproximationMax)



#objective
objective = myNetwork.Objective()
if (printResult):
    print 'objective including regularization terms (transformed space):\n  ' + str(objective)
	
#gradient
NumParameters = myNetwork.GetNumParameters()
gradient = nuto.DoubleFullMatrix(NumParameters,1)
myNetwork.Gradient(gradient)
if (printResult):
    print 'gradient:'
    gradient.Trans().Info()

#hessian
hessian = nuto.DoubleFullMatrix(NumParameters,NumParameters)
myNetwork.HessianFull(hessian)
if (printResult):
    print 'hessian:'
    hessian.Info()
