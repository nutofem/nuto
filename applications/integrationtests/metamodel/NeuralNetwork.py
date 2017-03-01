import nuto
import random
import math
import sys
import os
import numpy as np

# if set to true, the result will be generated (for later use in the test routine)
# otherwise, the current result will be compared to the stored result
createResult = False

# show the results on the screen
printResult = True

# system name and processor
system = sys.argv[1]+sys.argv[2]

# path in the original source directory and current filename at the and
pathToResultFiles = os.path.join(
        sys.argv[3], "results", system, os.path.basename(sys.argv[0]))

# remove the extension
fileExt = os.path.splitext(sys.argv[0])[1]
pathToResultFiles = pathToResultFiles.replace(fileExt, '')

# no error in file, modified, if error is detected
error = False

# Random seed for training data generation and initial weights
randomSeed = 1234567
random.seed(randomSeed)

# Generate Training data in the range min max
minCoordinate = -3
maxCoordinate = 3
meanNoise = 0
sigmaNoise = 5e-1
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
myNetwork.SetShowSteps(1)
# this is not very practical, but only used for the test file
myNetwork.SetMaxFunctionCalls(0)
myNetwork.SetVerboseLevel(verboseLevel)
myNetwork.SetMaxBayesianIterations(1)
myNetwork.UseDiagHessian()

myNetwork.SetSupportPoints(dimInput, dimOutput, SupportPointsInput, SupportPointsOutputWithNoise)

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

if (printResult):
    Parameters = myNetwork.GetParameters()
    print ('final set of parameters:\n', Parameters)

# calculate error for the training samples including standard deviation of
# approximation
SupportPointsApproximation = np.zeros((dimOutput, numSamples))
SupportPointsApproximationMin = np.zeros((dimOutput, numSamples))
SupportPointsApproximationMax = np.zeros((dimOutput, numSamples))
myNetwork.SolveConfidenceInterval(SupportPointsInput,
                                  SupportPointsApproximation,
                                  SupportPointsApproximationMin,
                                  SupportPointsApproximationMax)

objective = myNetwork.Objective()
if (printResult):
    print ('objective including regularization terms (transformed space):\n  ', objective)
if (createResult):
    f = open(pathToResultFiles+'Objective.txt', 'w')
    f.write('#Correct Objective\n')
    f.write(str(objective))
    f.close()
else:
    f = open(pathToResultFiles+'Objective.txt', 'r')
    f.readline()
    objectiveExact = float(f.readline())
    f.close()
    if (math.fabs(objective - objectiveExact) > 1e-8):
        print ('[' + system, sys.argv[0] + '] : objective is not correct.(' + str(math.fabs(objective - objectiveExact)) + ')')
        print ('objectiveExact including regularization terms (transformed space):\n', objectiveExact)
        error = True

# gradient
NumParameters = myNetwork.GetNumParameters()
gradient = np.zeros((NumParameters, 1))
myNetwork.Gradient(gradient)
gradient = gradient.squeeze()
if (printResult):
    print ('gradient:\n', gradient)
if (createResult):
    np.savetxt(pathToResultFiles+"Gradient.txt", gradient, header="#Correct gradient matrix")
else:
    gradientExact = np.loadtxt(pathToResultFiles + "Gradient.txt", skiprows=1)
    if (np.max(np.abs(gradientExact - gradient)) > 1e-8):
        print ('[' + system, sys.argv[0] + '] : gradient is not correct.(' + str(np.max(np.abs(gradientExact - gradient))) + ')')
        print ('gradientExact:\n', gradientExact)
        error = True

# hessian
hessian = np.zeros((NumParameters, NumParameters))
myNetwork.HessianFull(hessian)
if (printResult):
    print ('hessian:\n', hessian)
if (createResult):
    np.savetxt(pathToResultFiles+"Hessian.txt", hessian, header="#Correct hessian matrix")
else:
    hessianExact = np.loadtxt(pathToResultFiles+"Hessian.txt", skiprows=1)
    RelError = (hessianExact - hessian) / hessianExact
    if (np.max(np.abs(RelError)) > 1e-8):
        print ('[' + system, sys.argv[0] + '] : hessian is not correct.(' + str(np.max(np.abs(RelError))) + ')')
        print ('hessianExact:\n', hessianExact)
        print ('relative Error:\n', RelError)
        error = True

plotResult = False

if plotResult:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import numpy

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = [1, 2, 3, 4,  5, 6, 7,  8, 9, 10]
    y = [5, 6, 2, 3, 13, 4, 1,  2, 4,  8]
    z = [2, 3, 3, 3,  5, 7, 9, 11, 9, 10]

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    ax.scatter(SupportPointsInput[0, :], SupportPointsInput[1, :],
               SupportPointsOutputExact[1, :], c='r', marker='o')
    ax.scatter(SupportPointsInput[0, :], SupportPointsInput[1, :],
               SupportPointsOutputWithNoise[1, :], c='g', marker='o')
    ax.scatter(SupportPointsInput[0, :], SupportPointsInput[1, :],
               SupportPointsApproximation[1, :], c='b', marker='o')

    plt.show()

    programPause = raw_input("Press the <ENTER> key to continue...")


if (error):
    sys.exit(-1)
else:
    sys.exit(0)
