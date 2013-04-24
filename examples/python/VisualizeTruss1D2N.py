# -*- coding: utf-8 -*-
# $Id$
import nuto

try:
    #create structure
    myStructure = nuto.Structure(1)

    #create nodes
    myNode1 = myStructure.NodeCreate("displacements",nuto.DoubleFullVector((1,)))
    myStructure.NodeSetDisplacements(myNode1,nuto.DoubleFullVector((0.0,)))
    myNode2 = myStructure.NodeCreate("displacements",nuto.DoubleFullVector((6,)))
    myStructure.NodeSetDisplacements(myNode2,nuto.DoubleFullVector((0.3,)))
    myNode3 = myStructure.NodeCreate("displacements",nuto.DoubleFullVector((10,)))
    myStructure.NodeSetDisplacements(myNode3,nuto.DoubleFullVector((0.6,)))

    #create element
    myElement1 = myStructure.ElementCreate("Truss1D2N",nuto.IntFullVector((myNode1,myNode2)))
    myElement2 = myStructure.ElementCreate("Truss1D2N",nuto.IntFullVector((myNode3,myNode2)))

    #create constitutive law
    myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic")
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10)
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.1)

    #create section
    mySection1 = myStructure.SectionCreate ("TRUSS")
    myStructure.SectionSetArea(mySection1,0.01)

    myStructure.ElementSetIntegrationType(myElement1,"1D2NGauss2Ip","NOIPDATA")
    myStructure.ElementSetConstitutiveLaw(myElement1,myMatLin)
    myStructure.ElementSetSection(myElement1,mySection1)
    myStructure.ElementSetConstitutiveLaw(myElement2,myMatLin)
    myStructure.ElementSetSection(myElement2,mySection1)

    #visualize element
    myStructure.AddVisualizationComponentDisplacements()
    myStructure.AddVisualizationComponentEngineeringStrain()
    myStructure.AddVisualizationComponentEngineeringStress()
    myStructure.ExportVtkDataFile("Truss1D2N.vtk")
except nuto.Exception, e:
    print e.ErrorMessage()
