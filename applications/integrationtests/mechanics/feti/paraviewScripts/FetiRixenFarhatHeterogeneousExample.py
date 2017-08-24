# load the simple module
from paraview.simple import *

Disconnect()
Connect()
value = 'EngineeringStress'
for i in range(0,16):
    # open data file with corresponding reader
    reader = OpenDataFile("../FetiRixenFarhatHeterogeneousExampleResultDir_" + str(i) + "/Group9999_ElementsAll.pvd")
    reader.UpdatePipeline()

    # store time steps in vector
    timeSteps = reader.TimestepValues

    # create a new 'Calculator'
    calculator 	= Calculator()
    calculator.CoordinateResults = 1
    calculator.Function = 'coords+Displacements'

    # create a new 'Cell Data to Point Data'
    cellDatatoPointData1 = CellDatatoPointData(Input=calculator)

    # look at the last time step
    animationScene = GetAnimationScene()
    animationScene.AnimationTime = timeSteps[-1]

    # select representation
    displayProperties = GetDisplayProperties()
    displayProperties.Representation = 'Surface With Edges'

    view = GetActiveViewOrCreate('RenderView')

    # change background color to white
    view.Background = [1,1,1]

    # remove center axes visability
    view.CenterAxesVisibility = 0

    # remove orientation axes visability
    view.OrientationAxesVisibility = 0

    #set image size [width, height]
    #view.ViewSize = [ 1000, 2000 ]

    # get color transfer function/color map for 'EngineeringStress'
    colorTransferFunction = GetColorTransferFunction(value)
    # select a vector component or comment line to display magnitude
    #colorTransferFunction.VectorComponent = 0
    #colorTransferFunction.VectorMode = 'Component'
    # rescale color and/or opacity maps used to include current data range
    #colorTransferFunction.RescaleTransferFunction(-20,2)

    #displayProperties.LookupTable.VectorComponent = 6
    #displayProperties.RescaleTransferFunctionToDataRange(True)
    #displayProperties.LookupTable = MakeBlueToRedLT(-1, 1)

    # careful! EngineeringStress is now a point data due to the filter
    ColorBy(displayProperties, ('POINTS', value))
    displayProperties.RescaleTransferFunctionToDataRange(False)

    # show color bar/color legend
    displayProperties.SetScalarBarVisibility(view, True)
    scalarbar = GetScalarBar(displayProperties.LookupTable,view)
    scalarbar.Orientation ='Horizontal'
    # select position: bottom left is [0,0] top right is [1,1]
    scalarbar.Position = [0.5, 0.2]

    # select camera options
    camera = GetActiveCamera()
    camera.SetPosition(1,1,1)
    camera.SetFocalPoint(1,1,-1)
    camera.SetViewUp(0,1,0)
    ResetCamera()
    # important for 3D postprocessing
    camera.Roll(0)
    camera.Elevation(0)
    camera.Azimuth(0)

    Show()
    Render()
