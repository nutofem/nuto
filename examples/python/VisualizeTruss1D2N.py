# -*- coding: utf-8 -*-
import nuto
import numpy as np

try:
    # create structure
    myStructure = nuto.Structure(1)

    # create nodes
    myNode1 = myStructure.NodeCreateDOFs("displacements", np.array([1.0]))
    myStructure.NodeSetDisplacements(myNode1, np.array([0.0]))
    myNode2 = myStructure.NodeCreateDOFs("displacements", np.array([6.0]))
    myStructure.NodeSetDisplacements(myNode2, np.array([0.3]))
    myNode3 = myStructure.NodeCreateDOFs("displacements", np.array([10.0]))
    myStructure.NodeSetDisplacements(myNode3, np.array([0.6]))

    # create interpolation type
    myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D")
    myStructure.InterpolationTypeAdd(myInterpolationType, "coordinates", "equidistant1")
    myStructure.InterpolationTypeAdd(myInterpolationType, "displacements", "equidistant1")

    # create element
    myElement1 = myStructure.ElementCreate(myInterpolationType, [myNode1, myNode2])
    myElement2 = myStructure.ElementCreate(myInterpolationType, [myNode2, myNode3])

    # create constitutive law
    myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress")
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, "Youngs_Modulus", 10)
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, "Poissons_Ratio", 0.1)

    # create section
    mySection1 = myStructure.SectionCreate("TRUSS")
    myStructure.SectionSetArea(mySection1, 0.01)

    myStructure.ElementSetConstitutiveLaw(myElement1, myMatLin)
    myStructure.ElementSetSection(myElement1, mySection1)
    myStructure.ElementSetConstitutiveLaw(myElement2, myMatLin)
    myStructure.ElementSetSection(myElement2, mySection1)

    # visualize element
    visualizationGroup = myStructure.GroupCreate("Elements")
    myStructure.GroupAddElementsTotal(visualizationGroup)

    myStructure.AddVisualizationComponent(visualizationGroup, "Displacements")
    myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStrain")
    myStructure.AddVisualizationComponent(visualizationGroup, "EngineeringStress")
    myStructure.ExportVtkDataFileElements("Truss1D2N.vtk")

except nuto.Exception, e:
    print e.ErrorMessage()
