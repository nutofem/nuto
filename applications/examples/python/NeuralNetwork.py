import nuto
import random
import numpy as np

# Random seed for training data generation and initial weights
randomSeed = 1234567
random.seed(randomSeed)

# Generate Training data in the range min max
minCoordinate = -3
maxCoordinate = 3
meanNoise = 0
sigmaNoise = 1e-1
dimInput = 2
dimOutput = 2
numSamples = 10
numNeuronsHiddenLayer = 3
verboseLevel = 0


# some auxiliary routines
def RandomSamplesUniform(minCoordinate, maxCoordinate, dim, numSamples):
    sample = np.zeros((dim, numSamples))
    for count1 in range(0, dim):
        for count2 in range(0, numSamples):
            sample[count1, count2] = random.uniform(minCoordinate, maxCoordinate)
    return sample


def RandomSamplesGauss(mean, sigma, dim, numSamples):
    sample = np.zeros((dim, numSamples))
    for count1 in range(0, dim):
        for count2 in range(0, numSamples):
            sample[count1, count2] = random.gauss(mean, sigma)
    return sample


def ExactFunction(dimOutput, coordinates):
    _, numSamples = coordinates.shape
    sample = np.zeros((dimOutput, numSamples))
    for countdim in range(0, dimOutput):
        for count in range(0, numSamples):
            sample[countdim, count] = coordinates[(countdim+1) % dimInput, count]
    return sample

# create support points
SupportPointsInput = RandomSamplesUniform(minCoordinate, maxCoordinate, dimInput, numSamples)
SupportPointsOutputExact = ExactFunction(dimOutput, SupportPointsInput)
SupportPointsOutputNoise = RandomSamplesGauss(meanNoise, sigmaNoise, dimOutput, numSamples)
SupportPointsOutputWithNoise = SupportPointsOutputExact + SupportPointsOutputNoise

# create network with certain number of neurons in each hiddenlayer
myNetwork = nuto.NeuralNetwork([numNeuronsHiddenLayer])
myNetwork.InitRandomNumberGenerator(randomSeed)
myNetwork.SetBayesianTraining()
myNetwork.SetInitAlpha(1e-5)
myNetwork.SetAccuracyGradient(1e-5)
myNetwork.SetMinDeltaObjectiveBetweenRestarts(1e-4)
myNetwork.SetMinDeltaObjectiveBayesianIteration(1e-3)
myNetwork.SetShowSteps(100)
myNetwork.SetMaxFunctionCalls(1000000)
myNetwork.SetVerboseLevel(verboseLevel)
myNetwork.SetMaxBayesianIterations(1)
myNetwork.UseDiagHessian()

myNetwork.SetSupportPoints(dimInput, dimOutput, SupportPointsInput,
                           SupportPointsOutputWithNoise)

# set transformations for input/output
lowerBound = -1
upperBound = 1
myNetwork.AppendMinMaxTransformationInput(lowerBound, upperBound)
myNetwork.AppendMinMaxTransformationOutput(lowerBound, upperBound)
myNetwork.BuildTransformation()
# set transfer functions
myNetwork.SetTransferFunction(0, nuto.NeuralNetwork.TanSig)
myNetwork.SetTransferFunction(1, nuto.NeuralNetwork.PureLin)

# build/train the network
myNetwork.Build()

# calculate error for the training samples including standard deviation of
# approximation
SupportPointsApproximation = np.zeros((dimOutput, numSamples))
SupportPointsApproximationMin = np.zeros((dimOutput, numSamples))
SupportPointsApproximationMax = np.zeros((dimOutput, numSamples))
myNetwork.SolveConfidenceInterval(SupportPointsInput,
                                  SupportPointsApproximation,
                                  SupportPointsApproximationMin,
                                  SupportPointsApproximationMax)


# objective
print ('objective including regularization terms (transformed space):\n', myNetwork.Objective())
	
# gradient
NumParameters = myNetwork.GetNumParameters()
gradient = np.zeros((NumParameters, 1))
myNetwork.Gradient(gradient)
print ('gradient:\n', gradient.squeeze())

# hessian
hessian = np.zeros((NumParameters, NumParameters))
myNetwork.HessianFull(hessian)
print ('hessian:\n', hessian)
