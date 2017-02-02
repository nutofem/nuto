#to be run with pvpython (from paraview in its binary directory)
from paraview.simple import *
servermanager.LoadState("NonlocalDamagePlasticityModel.pvsm")
views = GetRenderViews()
textwidth=1920
textheight=1080
width=textwidth/3.-30
height=textheight/2.-60
pos_x=0
pos_y=0
for view in views:
    print pos_x , pos_y
    view.ViewSize = [width,height]
    window = view.GetRenderWindow()
    window.SetPosition(pos_x, pos_y)
    view.ViewPosition = [pos_x,pos_y]
    view.StillRender()

    if (pos_x+2*(width+10)<textwidth):
        pos_x+=width+10
    else:
        pos_x=0
	pos_y+=height+30
raw_input("Press ENTER to exit")




#from libvtkRenderingPython import *
#iren = vtkRenderWindowInteractor()
#iren.SetInteractorStyle(vtkInteractorStyleSwitch())
#iren.Initialize()
#iren.SetRenderWindow(view.GetRenderWindow())
#iren.Start()
