#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#define PRINTRESULT false

void CoefficientCheckLinearElasticTriangle(NuTo::Interpolation::eTypeOrder rTypeOrder)
{
    //create structure
    NuTo::Structure myStructure(2);
    //create nodes
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoordinates(2, 8);
    nodeCoordinates <<
            0, 10, 2, 8, 4, 8, 0, 10,
            0, 0, 2, 3, 7, 7, 10, 10;
    myStructure.NodesCreate(nodeCoordinates);

    int interpolationType = myStructure.InterpolationTypeCreate("Triangle2D");
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, rTypeOrder);

    //create element
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> nodeNumbers(3, 10);
    nodeNumbers <<
            0, 0, 0, 1, 2, 2, 3, 3, 4, 5,
            1, 2, 3, 7, 4, 3, 5, 7, 5, 7,
            3, 6, 2, 3, 6, 4, 4, 5, 6, 6;

    myStructure.ElementsCreate(interpolationType, nodeNumbers);
    myStructure.ElementTotalConvertToInterpolationType();


    //Calculate maximum independent sets for parallelization (openmp)
    myStructure.CalculateMaximumIndependentSets();

    //create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);

    //create section
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection, 1);

    //assign constitutive law
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection);

    // add a rather random constraint to get some dependent dofs
    NuTo::FullVector<double, Eigen::Dynamic> direction(2);
    direction(0) = 1;
    direction(1) = 0;
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    myStructure.NodeBuildGlobalDofs();

    bool isCorrectStiffnessStructure = myStructure.CheckCoefficientMatrix_0         (1.e-6, true);
    bool isCorrectStiffnessElements  = myStructure.ElementCheckCoefficientMatrix_0  (1.e-6);

    if (!isCorrectStiffnessStructure)
        throw NuTo::MechanicsException("Global stiffness matrix is incorrect.");

    if (!isCorrectStiffnessElements)
        throw NuTo::MechanicsException("Element stiffness matrices are incorrect.");
}

int main()
{
    try
    {
        CoefficientCheckLinearElasticTriangle(NuTo::Interpolation::EQUIDISTANT3);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return -1;
    }

    return 0;
}
