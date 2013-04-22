// $Id$

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
    // create structure
    NuTo::Structure myStructure(3);

    // create nodes
    NuTo::FullVector<double,Eigen::Dynamic> Coordinates(3);
    NuTo::FullVector<int,Eigen::Dynamic> Incidence(10);

    //create nodes
    Coordinates(0) = 0.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.0;
    Incidence(0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 1.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.0;
    Incidence(1) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 1.0;
    Coordinates(2) = 0.0;
    Incidence(2) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 1.0;
    Incidence(3) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.5;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.0;
    Incidence(4) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.5;
    Coordinates(1) = 0.5;
    Coordinates(2) = 0.0;
    Incidence(5) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.5;
    Coordinates(2) = 0.0;
    Incidence(6) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.5;
    Incidence(7) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.5;
    Coordinates(1) = 0.0;
    Coordinates(2) = 0.5;
    Incidence(8) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0) = 0.0;
    Coordinates(1) = 0.5;
    Coordinates(2) = 0.5;
    Incidence(9) = myStructure.NodeCreate("displacements",Coordinates);

	// create element
    int myElement1 = myStructure.ElementCreate("Tetrahedron10N",Incidence);

	//create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElastic");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,1);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0);

    // assign constitutive law
    myStructure.ElementSetConstitutiveLaw(myElement1,myMatLin);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> Ke;
    NuTo::FullVector<int,Eigen::Dynamic> rowIndex;
    NuTo::FullVector<int,Eigen::Dynamic> colIndex;
    myStructure.ElementStiffness(myElement1,Ke,rowIndex,colIndex);

	return 0;
}
