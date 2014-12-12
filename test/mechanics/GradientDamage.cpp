#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"

#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>




int main()
{
try {

    // define geometry
    int numElements = 7;
    double length = 50;
    double xWeakSpot = 25;
    double lWeakSpot = 5;

    double area = 10;
    double alpha = 0.10;
    double exponent = 4;


    // 1D structure
    NuTo::Structure myStructure(1);

    // create a damage law
    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(3);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;
    myDamageLaw(1) = 1.e-4;
    myDamageLaw(2) = 0.005;

    int myNumberConstitutiveLaw = myStructure.ConstitutiveLawCreate("GradientDamageEngineeringStress");
    myStructure.ConstitutiveLawSetDensity       (myNumberConstitutiveLaw, 1.0);
    myStructure.ConstitutiveLawSetYoungsModulus (myNumberConstitutiveLaw, 20000.);
    myStructure.ConstitutiveLawSetPoissonsRatio (myNumberConstitutiveLaw, 0.3);
    myStructure.ConstitutiveLawSetNonlocalRadius(myNumberConstitutiveLaw, 1.0);
    myStructure.ConstitutiveLawSetDamageLaw(myNumberConstitutiveLaw, myDamageLaw);

    // create sections
    int mySection    = myStructure.SectionCreate("Truss");
    myStructure.SectionSetArea             (mySection, area);
    myStructure.SectionSetDOF              (mySection, "displacements nonlocalEqStrain");
    myStructure.SectionSetInputConstitutive(mySection, "deformationGradient nonlocalEqStrain");

    // set function for area reduction
    NuTo::SectionTruss* secTruss = myStructure.SectionGetSectionPtr(mySection)->AsSectionTruss();
    double areaParameters[4];
    areaParameters[0] = xWeakSpot;
    areaParameters[1] = lWeakSpot;
    areaParameters[2] = alpha;
    areaParameters[3] = exponent;
    secTruss->SetAreaParameters(areaParameters);

    // create nodes
    int numNodes= numElements * 2 + 1;
    double lengthElement=length / numElements;

    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
    for(int node = 0; node < numNodes; node++)
    {
        nodeCoordinates(0) = node * lengthElement / 2.; // three nodes per element
        if (node % 2 == 0)
            myStructure.NodeCreate(node, "displacements nonlocalEqStrain", nodeCoordinates, 0);
        else
            myStructure.NodeCreate(node, "displacements", nodeCoordinates, 0);
    }

    int nodeLeft = 0;
    int nodeRight = numNodes - 1;

    // create elements

    for (int iElement = 0; iElement < numElements; iElement ++)
    {
        NuTo::FullVector<int, 3> nodeIds;
        nodeIds(0) = 2*iElement;
        nodeIds(1) = 2*iElement+1;
        nodeIds(2) = 2*iElement+2;
        myStructure.ElementCreate(iElement, "Truss1D3N", nodeIds,"ConstitutiveLawIpNonlocal","StaticData");
        myStructure.ElementSetIntegrationType(iElement,"1D2NGauss2Ip","StaticData");
        myStructure.ElementSetSection(iElement,mySection);
        myStructure.ElementSetConstitutiveLaw(iElement,myNumberConstitutiveLaw);
    }

    //     Constraints: set derivatives of nonlocal eq strain = 0
    int constraintEqLeft = myStructure.ConstraintLinearEquationCreate(nodeLeft, "nonlocalEqStrain", 1., 0.);
    myStructure.ConstraintLinearEquationAddTerm(constraintEqLeft, nodeLeft+2, "nonlocalEqStrain", -1);

    int constraintEqRight = myStructure.ConstraintLinearEquationCreate(nodeRight, "nonlocalEqStrain", 1., 0.);
    myStructure.ConstraintLinearEquationAddTerm(constraintEqRight, nodeRight-2, "nonlocalEqStrain", -1);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    bool globalStiffnessCorrect = myStructure.CheckCoefficientMatrix_0(1.e-6, true);
    bool elementStiffnessCorrect = myStructure.ElementCheckCoefficientMatrix_0(1.e-6);

    if (not globalStiffnessCorrect)
        throw NuTo::Exception("global stiffness matrix incorrect!");
    if (not elementStiffnessCorrect)
        throw NuTo::Exception("element stiffness matrices incorrect!");



} catch (NuTo::Exception& e) {
    std::cout << "Error executing GradientDamage "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    return -1;
}
    std::cout << "GradientDamage terminated normally."<< std::endl;
    return 0;
}
