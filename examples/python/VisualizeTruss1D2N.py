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
    myStructure.ConstitutiveLawCreate("myMatLin","LinearElastic")
    myStructure.ConstitutiveLawSetYoungsModulus("myMatLin",10)
    myStructure.ConstitutiveLawSetPoissonsRatio("myMatLin",0.1)

    #create section
    myStructure.SectionCreate("mySection1","1D")
    myStructure.SectionSetArea("mySection1",0.01)

    myStructure.ElementSetIntegrationType(myElement1,"1D2NGauss2Ip")
    myStructure.ElementSetConstitutiveLaw(myElement1,"myMatLin")
    myStructure.ElementSetSection(myElement1,"mySection1")
    myStructure.ElementSetConstitutiveLaw(myElement2,"myMatLin")
    myStructure.ElementSetSection(myElement2,"mySection1")

    #visualize element
    myStructure.ExportVtkDataFile("Truss1D2N.vtk","displacements engineering_strain engineering_stress")
except nuto.Exception, e:
    print e.ErrorMessage()