// $Id$

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
    // create structure
    NuTo::Structure myStructure(3);

    // create nodes
    NuTo::FullMatrix<double> Coordinates(3,1);
    NuTo::FullMatrix<int> Incidence(10,1);

    //create nodes
    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 0.0;
    Coordinates(2,0) = 0.0;
    Incidence(0,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 1.0;
    Coordinates(1,0) = 0.0;
    Coordinates(2,0) = 0.0;
    Incidence(1,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 1.0;
    Coordinates(2,0) = 0.0;
    Incidence(2,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 0.0;
    Coordinates(2,0) = 1.0;
    Incidence(3,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.5;
    Coordinates(1,0) = 0.0;
    Coordinates(2,0) = 0.0;
    Incidence(4,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.5;
    Coordinates(1,0) = 0.5;
    Coordinates(2,0) = 0.0;
    Incidence(5,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 0.5;
    Coordinates(2,0) = 0.0;
    Incidence(6,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 0.0;
    Coordinates(2,0) = 0.5;
    Incidence(7,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.5;
    Coordinates(1,0) = 0.0;
    Coordinates(2,0) = 0.5;
    Incidence(8,0) = myStructure.NodeCreate("displacements",Coordinates);

    Coordinates(0,0) = 0.0;
    Coordinates(1,0) = 0.5;
    Coordinates(2,0) = 0.5;
    Incidence(9,0) = myStructure.NodeCreate("displacements",Coordinates);

	// create element
    int myElement1 = myStructure.ElementCreate("Tetrahedron10N",Incidence);

	//create constitutive law
    myStructure.ConstitutiveLawCreate("myMatLin","LinearElastic");
    myStructure.ConstitutiveLawSetYoungsModulus("myMatLin",1);
    myStructure.ConstitutiveLawSetPoissonsRatio("myMatLin",0);

    // assign constitutive law
    myStructure.ElementSetConstitutiveLaw(myElement1,"myMatLin");

    NuTo::FullMatrix<double> Ke;
    NuTo::FullMatrix<int> rowIndex;
    NuTo::FullMatrix<int> colIndex;
    myStructure.ElementStiffness(myElement1,Ke,rowIndex,colIndex);

	return 0;
}
