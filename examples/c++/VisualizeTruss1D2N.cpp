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
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Coordinates(1,1);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Displacements(1,1);

    Coordinates(0,0) = 1;
    int myNode1 = myStructure.NodeCreate("displacements", Coordinates);
    Displacements(0,0) = 0;
    myStructure.NodeSetDisplacements(myNode1, Displacements);

    Coordinates(0,0) = 6;
    int myNode2 = myStructure.NodeCreate("displacements", Coordinates);
    Displacements(0,0) = 0.3;
    myStructure.NodeSetDisplacements(myNode2, Displacements);

    Coordinates(0,0) = 10;
    int myNode3 = myStructure.NodeCreate("displacements", Coordinates);
    Displacements(0,0) = 0.6;
    myStructure.NodeSetDisplacements(myNode3, Displacements);

    // create element
    NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> Incidences(2,1);

    Incidences(0,0) = myNode1;
    Incidences(1,0) = myNode2;
    int myElement1 = myStructure.ElementCreate("Truss1D2N", Incidences);

    Incidences(0,0) = myNode3;
    Incidences(1,0) = myNode2;
    int myElement2 = myStructure.ElementCreate("Truss1D2N", Incidences);

    // create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.1);

    // create section
    int mySection1 = myStructure.SectionCreate("1D");
    myStructure.SectionSetArea(mySection1,0.01);

    // assign material, section and integration type
    myStructure.ElementSetIntegrationType(myElement1,"1D2NGauss2Ip","NOIPDATA");
    myStructure.ElementSetConstitutiveLaw(myElement1,myMatLin);
    myStructure.ElementSetSection(myElement1,mySection1);
    myStructure.ElementSetConstitutiveLaw(myElement2,myMatLin);
    myStructure.ElementSetSection(myElement2,mySection1);

    // visualize element
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.ExportVtkDataFile("Truss1D2N.vtk");

    return 0;
}
