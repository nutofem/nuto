#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"

#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>




int main()
{
try {

    // define geometry
    int numElements = 3;
    double length = 50;

    double area = 10;
    double areaRed = area ;//* (1-0.10);

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
    int mySectionRed = myStructure.SectionCreate("Truss");
    myStructure.SectionSetArea             (mySection   , area);
    myStructure.SectionSetArea             (mySectionRed, areaRed);
    myStructure.SectionSetDOF              (mySection,    "displacements nonlocalEqStrain");
    myStructure.SectionSetDOF              (mySectionRed, "displacements nonlocalEqStrain");
    myStructure.SectionSetInputConstitutive(mySection,    "deformationGradient nonlocalEqStrain");
    myStructure.SectionSetInputConstitutive(mySectionRed, "deformationGradient nonlocalEqStrain");

    // create nodes
    int numNodes= numElements * 2 + 1;

    double lengthElement=length / numElements;

    // create nodes
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
    for(int node = 0; node < numNodes; node++)
    {
        nodeCoordinates(0) = node * lengthElement / 2.; // three nodes per element
        if (node % 2 == 0)
            myStructure.NodeCreate(node, "displacements nonlocalEqStrain", nodeCoordinates, 0);
        else
            myStructure.NodeCreate(node, "displacements", nodeCoordinates, 0);
    }

    myStructure.NodeBuildGlobalDofs();

    // create elements

    for (int element = 0; element < numElements; element ++)
    {
        NuTo::FullVector<int, 3> nodeIds;
        nodeIds(0) = 2*element;
        nodeIds(1) = 2*element+1;
        nodeIds(2) = 2*element+2;
        myStructure.ElementCreate(element, "Truss1D3N", nodeIds,"ConstitutiveLawIpNonlocal","StaticData");
        myStructure.ElementSetIntegrationType(element,"1D2NGauss2Ip","StaticData");
        if (element==numElements/2)
            myStructure.ElementSetSection(element,mySectionRed);
        else
            myStructure.ElementSetSection(element,mySection);
        myStructure.ElementSetConstitutiveLaw(element,myNumberConstitutiveLaw);
    }

    myStructure.CalculateMaximumIndependentSets();



    NuTo::SparseMatrixCSRVector2General<double> K_sparse;
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
    NuTo::FullVector<double,Eigen::Dynamic> activeDOFsCheck;
    NuTo::FullVector<double,Eigen::Dynamic> dependentDOFsCheck;
    myStructure.BuildGlobalCoefficientMatrix0(K_sparse, dispForceVector);

    myStructure.NodeExtractDofValues(0, activeDOFsCheck, dependentDOFsCheck);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> K(K_sparse);

    double epsilon = 1e-6;

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> K_num(K_sparse.GetNumRows(), K_sparse.GetNumColumns()), K_diff;
    NuTo::FullVector<double,Eigen::Dynamic> intForceVector1, intForceVector2;

    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector1);
    for (int count=0; count<activeDOFsCheck.GetNumRows(); count++)
    {
        activeDOFsCheck(count,0)+=epsilon;
        myStructure.NodeMergeActiveDofValues(0,activeDOFsCheck);
        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector2);
        K_num.SetColumn(count,(intForceVector2-intForceVector1)*(1./epsilon));
        activeDOFsCheck(count,0)-=epsilon;
    }
    myStructure.NodeMergeActiveDofValues(0,activeDOFsCheck);

    K_diff = K-K_num;
    double maxDiffK = abs(K_diff.Max());
    if (maxDiffK > 1e-3)
    {
        std::cout << "Something with the stiffness matrix, see values: " << std::endl;
        std::cout << "K:" << std::endl;
        K.Info(10,3, true);
        std::cout << "K_num:" << std::endl;
        K_num.Info(10,3, true);
        std::cout << "K - K_num:" << std::endl;
        K_diff.Info(10,3,true);
    } else {
        std::cout << "Stiffness matrix is kinda allright. Could be better, I guess. Nah.. it's allright." << std::endl;
        std::cout << "Max. difference: " << maxDiffK << std::endl;

    }



} catch (NuTo::Exception& e) {
    std::cout << "Error executing GradientDamage "<< std::endl;
    std::cout << e.ErrorMessage() << std::endl;
    return -1;
}
    std::cout << "GradientDamage terminated normally."<< std::endl;
    return 0;
}
