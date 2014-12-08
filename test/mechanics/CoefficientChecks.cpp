#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#define PRINTRESULT false
#define SUBDIVIDE_ELEMENTS false

int main()
{
    try
    {
        //create structure
        NuTo::Structure myStructure(2);

        //create nodes
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoordinates(2,8);
        nodeCoordinates <<
                0, 10, 2, 8, 4, 8,  0, 10,
                0,  0, 2, 3, 7, 7, 10, 10;

        myStructure.NodesCreate("displacements",nodeCoordinates);

        //create element
        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> nodeNumbers(3,10);
        nodeNumbers <<
                0, 0, 0, 1, 2, 2, 3, 3, 4, 5,
                1, 2, 3, 7, 4, 3, 5, 7, 5, 7,
                3, 6, 2, 3, 6, 4, 4, 5, 6, 6;

        myStructure.ElementsCreate("PLANE2D3N", nodeNumbers);

        if (SUBDIVIDE_ELEMENTS)
        {
            //convert 3N into 10N elements
            int elementGroup = myStructure.GroupCreate("Elements");
            myStructure.GroupAddElementFromType(elementGroup, "PLANE2D3N");
            myStructure.ElementConvertPlane2D3N(elementGroup,"PLANE2D10N",1e-5,1.);
        }

        //Calculate maximum independent sets for parallelization (openmp)
        myStructure.CalculateMaximumIndependentSets();

        //create constitutive law
        int myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
        myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,10);
        myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,0.2);

        //create section
        int mySection = myStructure.SectionCreate("Plane_Strain");
        myStructure.SectionSetThickness(mySection,1);

        //assign constitutive law
        myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
        myStructure.ElementTotalSetSection(mySection);

        // add a rather random constraint to get some dependent dofs
        NuTo::FullVector<double,Eigen::Dynamic> direction(2);
        direction(0) = 1;
        direction(1) = 0;

        myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);

        myStructure.NodeBuildGlobalDofs();

        bool isCorrectStiffnessStructure = myStructure.CheckCoefficientMatrix_0(1.e-6, true);

        bool isCorrectStiffnessElements = myStructure.ElementCheckCoefficientMatrix_0(1.e-6);

        if (not isCorrectStiffnessStructure)
            throw NuTo::MechanicsException("Global stiffness matrix is incorrect.");
        if (not isCorrectStiffnessElements)
            throw NuTo::MechanicsException("Element stiffness matrices are incorrect.");



    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return -1;
    }

    return 0;
}
