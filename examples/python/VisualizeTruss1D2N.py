# -*- coding: utf-8 -*-
# $Id$
import nuto

try:
    #create structure
    myStructure = nuto.Structure(1)

    #create nodes
    myNode1 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(1,1,(1,)))
    myStructure.NodeSetDisplacements(myNode1,nuto.DoubleFullMatrix(1,1,(0.0,)))
    myNode2 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(1,1,(6,)))
    myStructure.NodeSetDisplacements(myNode2,nuto.DoubleFullMatrix(1,1,(0.3,)))
    myNode3 = myStructure.NodeCreate("displacements",nuto.DoubleFullMatrix(1,1,(10,)))
    myStructure.NodeSetDisplacements(myNode3,nuto.DoubleFullMatrix(1,1,(0.6,)))

    #create element
    myElement1 = myStructure.ElementCreate("Truss1D2N",nuto.IntFullMatrix(2,1,(myNode1,myNode2)))
    myElement2 = myStructure.ElementCreate("Truss1D2N",nuto.IntFullMatrix(2,1,(myNode3,myNode2)))

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
