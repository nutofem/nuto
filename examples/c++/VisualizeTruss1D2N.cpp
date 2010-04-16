// $Id$

#include "nuto/math/MathException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
    // create structure
    NuTo::Structure myStructure(1);

    // create nodes
    NuTo::FullMatrix<double> Coordinates(1,1);
    NuTo::FullMatrix<double> Displacements(1,1);

    Coordinates(0,0) = 1;
    unsigned int myNode1 = myStructure.NodeCreate("displacements", Coordinates);
    Displacements(0,0) = 0;
    myStructure.NodeSetDisplacements(myNode1, Displacements);

    Coordinates(0,0) = 6;
    unsigned int myNode2 = myStructure.NodeCreate("displacements", Coordinates);
    Displacements(0,0) = 0.3;
    myStructure.NodeSetDisplacements(myNode2, Displacements);

    Coordinates(0,0) = 10;
    unsigned int myNode3 = myStructure.NodeCreate("displacements", Coordinates);
    Displacements(0,0) = 0.6;
    myStructure.NodeSetDisplacements(myNode3, Displacements);

    // create element
    NuTo::FullMatrix<int> Incidences(2,1);

    Incidences(0,0) = myNode1;
    Incidences(1,0) = myNode2;
    unsigned int myElement1 = myStructure.ElementCreate("Truss1D2N", Incidences);

    Incidences(0,0) = myNode3;
    Incidences(1,0) = myNode2;
    unsigned int myElement2 = myStructure.ElementCreate("Truss1D2N", Incidences);

    // create constitutive law
    myStructure.ConstitutiveLawCreate("myMatLin","LinearElastic");
    myStructure.ConstitutiveLawSetYoungsModulus("myMatLin",10);
    myStructure.ConstitutiveLawSetPoissonsRatio("myMatLin",0.1);

    // create section
    myStructure.SectionCreate("mySection1","1D");
    myStructure.SectionSetArea("mySection1",0.01);

    // assign material, section and integration type
    myStructure.ElementSetIntegrationType(myElement1,"1D2NGauss2Ip");
    myStructure.ElementSetConstitutiveLaw(myElement1,"myMatLin");
    myStructure.ElementSetSection(myElement1,"mySection1");
    myStructure.ElementSetConstitutiveLaw(myElement2,"myMatLin");
    myStructure.ElementSetSection(myElement2,"mySection1");

    // visualize element
    myStructure.ExportVtkDataFile("Truss1D2N.vtk","displacements engineering_strain engineering_stress");

    return 0;
}
