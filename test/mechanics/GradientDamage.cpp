#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"

//#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>



bool CheckDamageLawsDerivatives(NuTo::GradientDamageEngineeringStress rConstitutiveLaw)
{
    double epsilon = 1.e-8;
    double E = rConstitutiveLaw.GetYoungsModulus();
    double e0 = rConstitutiveLaw.GetTensileStrength() / E;
    double step = e0/5;
    for (int i = 1; i < 100; ++i)
    {
        double kappa = i*step+epsilon;
//        kappa = i*step;
        double sigma1 = (1 - rConstitutiveLaw.CalculateDamage(kappa))*E*kappa;
        double sigma2 = (1 - rConstitutiveLaw.CalculateDamage(kappa+epsilon))*E*(kappa+epsilon);

        double DsigmaDkappa = -rConstitutiveLaw.CalculateDerivativeDamage(kappa)*E*kappa
                + (1-rConstitutiveLaw.CalculateDamage(kappa))*E;
        double DsigmaDkappa_CDF = (sigma2-sigma1)/epsilon;

        double differenceSigma = DsigmaDkappa - DsigmaDkappa_CDF;

#ifdef PRINTRESULT
            std::cout << "kappa:" <<  kappa << " | differenceSigma: " << differenceSigma << std::endl;
            std::cout << "Dsigma:" << DsigmaDkappa << " | Dsigma_CDF: " << DsigmaDkappa_CDF << std::endl;
            std::cout << "sigma1:" << sigma1 << " | sigma2: " << sigma2 << std::endl << std::endl;
#endif

        if (std::abs(differenceSigma)>1.e-3)
            return false;
    }

    return true;
}

void CheckDamageLaws()
{
    NuTo::GradientDamageEngineeringStress myConstitutiveLaw;


    myConstitutiveLaw.SetDensity        (1.0);
    myConstitutiveLaw.SetYoungsModulus  (30000);
    myConstitutiveLaw.SetPoissonsRatio  (0.3);
    myConstitutiveLaw.SetNonlocalRadius (1.0);
    myConstitutiveLaw.SetTensileStrength(4.);
    myConstitutiveLaw.SetFractureEnergy (0.21);

    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

    // create a damage law
    NuTo::FullVector<double, Eigen::Dynamic> myDamageLawNoSoftening(1);
    myDamageLawNoSoftening(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_NO_SOFTENING;


    NuTo::FullVector<double, Eigen::Dynamic> myDamageLawLinear(1);
    myDamageLawLinear(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_LINEAR_SOFTENING;

    NuTo::FullVector<double, Eigen::Dynamic> myDamageLawExponential(1);
    myDamageLawExponential(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

    NuTo::FullVector<double, Eigen::Dynamic> myDamageLawHermite(1);
    myDamageLawHermite(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_CUBIC_HERMITE;

//    NuTo::FullVector<double, Eigen::Dynamic> myDamageLawExponentialSmooth(4);
//    myDamageLawExponentialSmooth(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH;
//    myDamageLawExponentialSmooth(1) = 0;
//    myDamageLawExponentialSmooth(2) = 10;
//    myDamageLawExponentialSmooth(3) = 20;


    myConstitutiveLaw.SetDamageLaw(myDamageLawNoSoftening);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
            throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_NO_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetDamageLaw(myDamageLawLinear);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
            throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_LINEAR_SOFTENING: wrong damage derivatives");


    myConstitutiveLaw.SetDamageLaw(myDamageLawExponential);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
            throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_EXPONENTIAL_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetDamageLaw(myDamageLawHermite);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
            throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_CUBIC_HERMITE: wrong damage derivatives");

//    myConstitutiveLaw.SetDamageLaw(myDamageLawExponentialSmooth);
//    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
//            throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH: wrong damage derivatives");

}


int main()
{
//try {

    CheckDamageLaws();

    {
        // USE TRUSS1D3N ELEMENTS

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

        NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(4);
        myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_SMOOTH;
        myDamageLaw(1) = 5;
        myDamageLaw(2) = 10;
        myDamageLaw(3) = 20;
        // create a damage law
        int myNumberConstitutiveLaw = myStructure.ConstitutiveLawCreate("GradientDamageEngineeringStress");
        myStructure.ConstitutiveLawSetDensity        (myNumberConstitutiveLaw, 1.0);
        myStructure.ConstitutiveLawSetYoungsModulus  (myNumberConstitutiveLaw, 30000);
        myStructure.ConstitutiveLawSetPoissonsRatio  (myNumberConstitutiveLaw, 0.3);
        myStructure.ConstitutiveLawSetNonlocalRadius (myNumberConstitutiveLaw, 1.0);
        myStructure.ConstitutiveLawSetTensileStrength(myNumberConstitutiveLaw, 4.);
        myStructure.ConstitutiveLawSetFractureEnergy (myNumberConstitutiveLaw, 0.21);
        myStructure.ConstitutiveLawSetDamageLaw      (myNumberConstitutiveLaw, myDamageLaw);

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

    }




    {

        // USE TRUSS1D4N ELEMENTS

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

        NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
        myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;
        // create a damage law
        int myNumberConstitutiveLaw = myStructure.ConstitutiveLawCreate("GradientDamageEngineeringStress");
        myStructure.ConstitutiveLawSetDensity        (myNumberConstitutiveLaw, 1.0);
        myStructure.ConstitutiveLawSetYoungsModulus  (myNumberConstitutiveLaw, 30000);
        myStructure.ConstitutiveLawSetPoissonsRatio  (myNumberConstitutiveLaw, 0.3);
        myStructure.ConstitutiveLawSetNonlocalRadius (myNumberConstitutiveLaw, 1.0);
        myStructure.ConstitutiveLawSetTensileStrength(myNumberConstitutiveLaw, 4.);
        myStructure.ConstitutiveLawSetFractureEnergy (myNumberConstitutiveLaw, 0.21);
        myStructure.ConstitutiveLawSetDamageLaw      (myNumberConstitutiveLaw, myDamageLaw);

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
        int numNodes = numElements + 1;
        double lengthElement=length / numElements;

        NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
        for(int node = 0; node < numNodes; node++)
        {
            nodeCoordinates(0) = node * lengthElement;
            myStructure.NodeCreate(node, "displacements nonlocalEqStrain", nodeCoordinates, 0);
        }

        int nodeLeft = 0;
        int nodeRight = numNodes - 1;


        int elementGroup = myStructure.GroupCreate("Elements");
        // create elements
        for (int iElement = 0; iElement < numElements; iElement ++)
        {
            NuTo::FullVector<int, 2> nodeIds;
            nodeIds(0) = iElement;
            nodeIds(1) = iElement+1;
            myStructure.ElementCreate(iElement, "Truss1D2N", nodeIds,"ConstitutiveLawIp","StaticData");
            myStructure.ElementSetIntegrationType(iElement,"1D2NGauss3Ip","StaticData");
            myStructure.ElementSetSection(iElement,mySection);
            myStructure.ElementSetConstitutiveLaw(iElement,myNumberConstitutiveLaw);
            myStructure.GroupAddElement(elementGroup, iElement);
        }

        myStructure.ElementConvertTruss1D2NToTruss1D4NDisp3NX(elementGroup, "nonlocalEqStrain");


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

    }



//} catch (NuTo::Exception& e) {
//    std::cout << "Error executing GradientDamage "<< std::endl;
//    std::cout << e.ErrorMessage() << std::endl;
//    return -1;
//}
    std::cout << "GradientDamage terminated normally."<< std::endl;
    return 0;
}
