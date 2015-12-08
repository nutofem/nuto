#include <boost/filesystem.hpp>

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"

#include "nuto/mechanics/elements/BoundaryElement1D.h"
#include "nuto/mechanics/elements/BoundaryElement2D.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"

//#define PRINTRESULT
#include <eigen3/Eigen/Eigenvalues>

bool CheckDamageLawsDerivatives(NuTo::GradientDamageEngineeringStress rConstitutiveLaw)
{
    double epsilon = 1.e-8;
    double E = rConstitutiveLaw.GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS);
    double e0 = rConstitutiveLaw.GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH) / E;
    double step = e0 / 5;
    for (int i = 1; i < 100; ++i)
    {
        double kappa = i * step + epsilon;
//        kappa = i*step;
        double sigma1 = (1 - rConstitutiveLaw.CalculateDamage(kappa)) * E * kappa;
        double sigma2 = (1 - rConstitutiveLaw.CalculateDamage(kappa + epsilon)) * E * (kappa + epsilon);

        double DsigmaDkappa = -rConstitutiveLaw.CalculateDerivativeDamage(kappa) * E * kappa + (1 - rConstitutiveLaw.CalculateDamage(kappa)) * E;
        double DsigmaDkappa_CDF = (sigma2 - sigma1) / epsilon;

        double differenceSigma = DsigmaDkappa - DsigmaDkappa_CDF;

#ifdef PRINTRESULT
        std::cout << "kappa:" << kappa << " | differenceSigma: " << differenceSigma << std::endl;
        std::cout << "Dsigma:" << DsigmaDkappa << " | Dsigma_CDF: " << DsigmaDkappa_CDF << std::endl;
        std::cout << "sigma1:" << sigma1 << " | sigma2: " << sigma2 << std::endl << std::endl;
#endif

        if (std::abs(differenceSigma) > 1.e-3)
            return false;
    }

    return true;
}

void CheckDamageLaws()
{
    NuTo::GradientDamageEngineeringStress myConstitutiveLaw;

    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY,1.0);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,30000);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.3);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,1.0);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,4.);
    myConstitutiveLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,0.21);

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

    myConstitutiveLaw.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW,myDamageLawNoSoftening);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_NO_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW,myDamageLawLinear);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_LINEAR_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW,myDamageLawExponential);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_EXPONENTIAL_SOFTENING: wrong damage derivatives");

    myConstitutiveLaw.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW,myDamageLawHermite);
    if (not CheckDamageLawsDerivatives(myConstitutiveLaw))
        throw NuTo::MechanicsException("DamageLaw::ISOTROPIC_CUBIC_HERMITE: wrong damage derivatives");

}

void CheckLocalEqStrainDerivatives()
{
    NuTo::GradientDamageEngineeringStress law;
    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

    // create a damage law

    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY,1.0);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,30000);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.3);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,1.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,4.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH,40.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,0.21);
    law.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW,myDamageLaw);

    NuTo::EngineeringStrain2D strain;
    NuTo::LocalEqStrain localEqStrain0, localEqStrain1;
    NuTo::ConstitutiveTangentLocal<3, 1> tangent, dummy, tangent_CDF;

    std::vector<NuTo::FullVector<double, 3>> strainCases;

    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 1, 0, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, 1, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, 0, 1 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 1, 1, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, 1, 1 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 1, 0, 1 }));

    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 3, 0, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, 3, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, 0, 3 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 3, 2, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, 2, 3 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 2, 0, 3 }));

    strainCases.push_back(NuTo::FullVector<double, 3>(
    { -3, 0, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, -3, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, 0, -3 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { -3, 2, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 0, -2, 3 }));
    strainCases.push_back(NuTo::FullVector<double, 3>(
    { 2, 0, -3 }));

//    strainCases.push_back(NuTo::FullVector<double, 3>({ 0, 0, 0 }));

    double delta = 1.e-8;

    // check derivatives for plane stress and plane strain
    for (int iSectionType = 0; iSectionType < 2; ++iSectionType)
    {
        for (unsigned int iCase = 0; iCase < strainCases.size(); ++iCase)
        {
            strain[0] = strainCases[iCase][0];
            strain[1] = strainCases[iCase][1];
            strain[2] = strainCases[iCase][2];

            // calculate derivative numerically

            if (iSectionType == 0)
            {
                law.CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStrain(strain, localEqStrain0, tangent);
            } else
            {
                law.CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStress(strain, localEqStrain0, tangent);
            }

            for (int i = 0; i < 3; ++i)
            {
                strain[i] += delta;

                if (iSectionType == 0)
                {
                    law.CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStrain(strain, localEqStrain1, dummy);
                } else
                {
                    law.CalculateLocalEqStrainAndDerivativeModifiedMises2DPlaneStress(strain, localEqStrain1, dummy);
                }

                tangent_CDF[i] = (localEqStrain1[0] - localEqStrain0[0]) / delta;

                strain[i] -= delta;
            }

            if ((tangent - tangent_CDF).cwiseAbs().maxCoeff() > 1.e-6)
            {
                std::cout << "strain:" << std::endl;
                strain.Info(15, 5, true);

                std::cout << "tangent algo: ";
                tangent.Trans().Info(15, 5, true);

                std::cout << "tangent cdf : ";
                tangent_CDF.Trans().Info(15, 5, true);
                throw NuTo::MechanicsException("[CheckLocalEqStrainDerivatives] wrong derivatives!");

            }

        }
    }
}

void CheckLocalEqStrainDerivatives3D()
{
    NuTo::GradientDamageEngineeringStress law;
    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

    // create a damage law

    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY,1.0);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,30000);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0.3);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,1.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,4.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH,40.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,0.21);
    law.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW,myDamageLaw);

    NuTo::EngineeringStrain3D strain;
    NuTo::LocalEqStrain localEqStrain0, localEqStrain1;
    NuTo::ConstitutiveTangentLocal<6, 1> tangent, dummy, tangent_CDF;

    std::vector<NuTo::FullVector<double, 6>> strainCases;

    strainCases.push_back(NuTo::FullVector<double, 6>(
    { 1, 0, 0, 0, 0, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 6>(
    { 0, 1, 0, 0, 0, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 6>(
    { 0, 0, 1, 0, 0, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 6>(
    { 0, 0, 0, 1, 0, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 6>(
    { 0, 0, 0, 0, 1, 0 }));
    strainCases.push_back(NuTo::FullVector<double, 6>(
    { 0, 0, 0, 0, 0, 1 }));

    double delta = 1.e-8;

    // check derivatives for plane stress and plane strain
    for (unsigned int iCase = 0; iCase < strainCases.size(); ++iCase)
    {
        strain[0] = strainCases[iCase][0];
        strain[1] = strainCases[iCase][1];
        strain[2] = strainCases[iCase][2];
        strain[3] = strainCases[iCase][3];
        strain[4] = strainCases[iCase][4];
        strain[5] = strainCases[iCase][5];

        law.CalculateLocalEqStrainAndDerivativeModifiedMises3D(strain, localEqStrain0, tangent);
        // calculate derivative numerically
        for (int i = 0; i < 6; ++i)
        {
            strain[i] += delta;

            law.CalculateLocalEqStrainAndDerivativeModifiedMises3D(strain, localEqStrain1, dummy);

            tangent_CDF[i] = (localEqStrain1[0] - localEqStrain0[0]) / delta;

            strain[i] -= delta;
        }

        if ((tangent - tangent_CDF).cwiseAbs().maxCoeff() > 1.e-7)
        {
            std::cout << "strain:" << std::endl;
            strain.Info(15, 5, true);

            std::cout << "tangent algo: ";
            tangent.Trans().Info(15, 5, true);

            std::cout << "tangent cdf : ";
            tangent_CDF.Trans().Info(15, 5, true);
            throw NuTo::MechanicsException("[CheckLocalEqStrainDerivatives] wrong derivatives!");

        }

    }
}

int SetConstitutiveLaw(NuTo::Structure& rStructure)
{
    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

    // create a damage law
    int lawId = rStructure.ConstitutiveLawCreate("GradientDamageEngineeringStress");
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::DENSITY, 1.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 1.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.21);
    rStructure.ConstitutiveLawSetDamageLaw(lawId, myDamageLaw);

//    int myNumberConstitutiveLaw = rStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
//    rStructure.ConstitutiveLawSetDensity(myNumberConstitutiveLaw, 1.0);
//    rStructure.ConstitutiveLawSetYoungsModulus(myNumberConstitutiveLaw, 30000);
//    rStructure.ConstitutiveLawSetPoissonsRatio(myNumberConstitutiveLaw, 0.0);

    std::cout << "lawId " << lawId << std::endl;
    return lawId;
}

void GradientDamage1D()
{
    // define geometry
    const int numElements = 5;
    const double length = 100;
    const double xWeakSpot = 25;
    const double lWeakSpot = 5;

    const double area = 10;
    const double alpha = 0.10;
    const double exponent = 4;

    const double BoundaryDisplacement = 1;

    // 1D structure
    NuTo::Structure myStructure(1);
    myStructure.SetVerboseLevel(10);

    int myNumberConstitutiveLaw = SetConstitutiveLaw(myStructure);

    // create sections
    int mySection = myStructure.SectionCreate("Truss");
    myStructure.SectionSetArea(mySection, area);

    // set function for area reduction
    NuTo::SectionTruss* secTruss = myStructure.SectionGetSectionPtr(mySection)->AsSectionTruss();
    double areaParameters[4];
    areaParameters[0] = xWeakSpot;
    areaParameters[1] = lWeakSpot;
    areaParameters[2] = alpha;
    areaParameters[3] = exponent;
    secTruss->SetAreaParameters(areaParameters);

    // create nodes
    int numNodes = numElements + 1; // nodes for nonlocal strain/coordinates
    double lengthElement = length / numElements;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(1);
    for (int node = 0; node < numNodes; node++)
    {
        nodeCoordinates(0) = node * lengthElement; // two nodes per element
        myStructure.NodeCreate(node, nodeCoordinates);
    }

    int interpolationType = myStructure.InterpolationTypeCreate("TRUSS1D");
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);

    // create elements
    NuTo::FullVector<int, Eigen::Dynamic> elementIncidence(2);
    for (int iElement = 0; iElement < numElements; iElement++)
    {
        elementIncidence(0) = iElement;
        elementIncidence(1) = iElement + 1;
        myStructure.ElementCreate(iElement, interpolationType, elementIncidence, "ConstitutiveLawIp", "StaticData");
    }

    int groupElements = myStructure.GroupCreate("ELEMENTS");
    myStructure.GroupAddElementFromType(groupElements, interpolationType);
    myStructure.ElementConvertToInterpolationType(groupElements);

    myStructure.InterpolationTypeSetIntegrationType(interpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::STATICDATA);
    myStructure.ElementTotalSetConstitutiveLaw(myNumberConstitutiveLaw);
    myStructure.ElementTotalSetSection(mySection);

    // set boundary conditions and loads
    NuTo::FullVector<double, Eigen::Dynamic> direction(1);
    direction(0) = 1;

    int nodeGroupBoundary = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0, -1.e-6, 1.e-6);
    myStructure.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0, length - 1.e-6, length + 1.e-6);

    int nodeGroupLeft = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNodeCoordinateRange(nodeGroupLeft, 0, -1.e-6, 1.e-6);
    int nodeIndexLeft = myStructure.GroupGetMemberIds(nodeGroupLeft).GetValue(0);

    int nodeGroupRight = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNodeCoordinateRange(nodeGroupRight, 0, length - 1.e-6, length + 1.e-6);
    int nodeIndexRight = myStructure.GroupGetMemberIds(nodeGroupRight).GetValue(0);

    int elemGroupBoundary = myStructure.GroupCreate("ELEMENTS");
    myStructure.GroupAddElementsFromNodes(elemGroupBoundary, nodeGroupBoundary, false);

    assert(myStructure.GroupGetNumMembers(nodeGroupBoundary) == 2);
    assert(myStructure.GroupGetNumMembers(elemGroupBoundary) == 2);

    NuTo::BoundaryType::eType boundaryConditionType = NuTo::BoundaryType::ROBIN_INHOMOGENEOUS;
    int gBoundaryElements = myStructure.BoundaryElementsCreate(elemGroupBoundary, nodeGroupBoundary);
    auto boundaryElementIds = myStructure.GroupGetMemberIds(gBoundaryElements);
    for (int iBoundaryElement = 0; iBoundaryElement < boundaryElementIds.GetNumRows(); ++iBoundaryElement)
    {
        NuTo::BoundaryElement1D* boundaryElement = myStructure.ElementGetElementPtr(boundaryElementIds.GetValue(iBoundaryElement))->AsBoundaryElement1D();
        boundaryElement->SetBoundaryConditionType(boundaryConditionType);
    }

    // Constraints
    myStructure.ConstraintLinearSetDisplacementNode(nodeIndexLeft, direction, 0.0);

    // displacement controlled load
    myStructure.ConstraintLinearSetDisplacementNode(nodeIndexRight, direction, BoundaryDisplacement);

    // build global dof numbering
    myStructure.NodeBuildGlobalDofs();

//    myStructure.Info();
    myStructure.CalculateMaximumIndependentSets();

    // apply linear displacement state and some nonlocal eq strains
    int allNodes = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNodeCoordinateRange(allNodes, 0, -1, length + 1);
    auto nodeIds = myStructure.GroupGetMemberIds(allNodes);
    for (int i = 0; i < nodeIds.GetNumRows(); ++i)
    {
        NuTo::NodeBase* node = myStructure.NodeGetNodePtr(nodeIds.GetValue(i));
        auto disps = node->GetCoordinates1D() / length * BoundaryDisplacement;
        node->SetDisplacements1D(disps);

        if (node->GetNumNonlocalEqStrain() > 0)
            node->SetNonlocalEqStrain(disps.at(0, 0) / 10);
    }

    myStructure.Info();
    std::cout << "-----------------------------------------------------------------" << std::endl;

    bool globalStiffnessCorrect = myStructure.CheckCoefficientMatrix_0(1.e-8, true);
    bool elementStiffnessCorrect = myStructure.ElementCheckCoefficientMatrix_0(1.e-8);

    if (not globalStiffnessCorrect)
        throw NuTo::Exception("[GradientDamage1D] global stiffness matrix incorrect!");
    if (not elementStiffnessCorrect)
        throw NuTo::Exception("[GradientDamage1D] element stiffness matrices incorrect!");

}

void GradientDamage2D()
{
    // define geometry

    double lX = 20.;
    double lY = 5.;

    double lZ = 2.;

    int numElementsX = 4;
    int numElementsY = 1;

    const double BoundaryDisplacement = 1;

    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);
    myStructure.SetVerboseLevel(10);

    //create nodes
    int numNodesX = numElementsX + 1;
    int numNodesY = numElementsY + 1;
    double deltaX = lX / (numElementsX);
    double deltaY = lY / (numElementsY);

    int nodeNum = 0;
    for (int countY = 0; countY < numNodesY; countY++)
    {
        for (int countX = 0; countX < numNodesX; countX++)
        {
            NuTo::FullVector<double, Eigen::Dynamic> coordinates(2);
            coordinates(0) = countX * deltaX;
            coordinates(1) = countY * deltaY;
            myStructure.NodeCreate(nodeNum, coordinates);
            nodeNum++;
        }
    }

    int myInterpolationType = myStructure.InterpolationTypeCreate("Quad2D");
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);

    //create elements
    for (int countY = 0; countY < numElementsY; countY++)
    {
        for (int countX = 0; countX < numElementsX; countX++)
        {
            NuTo::FullVector<int, Eigen::Dynamic> nodes(4);
            nodes(0) = countX + countY * numNodesX;
            nodes(1) = countX + 1 + countY * numNodesX;
            nodes(2) = countX + 1 + (countY + 1) * numNodesX;
            nodes(3) = countX + (countY + 1) * numNodesX;
            myStructure.ElementCreate(myInterpolationType, nodes, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);
        }
    }

    std::cout << myStructure.ElementGetElementPtr(0)->GetNumIntegrationPoints() << std::endl;
//    return ;

    int myConstitutiveLaw = SetConstitutiveLaw(myStructure);
    int mySection = myStructure.SectionCreate("Plane_Strain");
    myStructure.SectionSetThickness(mySection, lZ);

    myStructure.ElementTotalConvertToInterpolationType();
    myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D4NGauss4Ip, NuTo::IpData::STATICDATA);
    myStructure.ElementTotalSetConstitutiveLaw(myConstitutiveLaw);
    myStructure.ElementTotalSetSection(mySection);

    int nodeGroupBoundary = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0, -1.e-6, 1.e-6);
    myStructure.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0, lX - 1.e-6, lX + 1.e-6);

    int elemGroupBoundary = myStructure.GroupCreate("ELEMENTS");
    myStructure.GroupAddElementsFromNodes(elemGroupBoundary, nodeGroupBoundary, false);

    assert(myStructure.GroupGetNumMembers(nodeGroupBoundary) >= 4);
    assert(myStructure.GroupGetNumMembers(elemGroupBoundary) == 2 * numElementsY);

    NuTo::BoundaryType::eType boundaryConditionType = NuTo::BoundaryType::ROBIN_INHOMOGENEOUS;
    int gBoundaryElements = myStructure.BoundaryElementsCreate(elemGroupBoundary, nodeGroupBoundary);
    auto boundaryElementIds = myStructure.GroupGetMemberIds(gBoundaryElements);
    for (int iBoundaryElement = 0; iBoundaryElement < boundaryElementIds.GetNumRows(); ++iBoundaryElement)
    {
        NuTo::BoundaryElement2D* boundaryElement = myStructure.ElementGetElementPtr(boundaryElementIds.GetValue(iBoundaryElement))->AsBoundaryElement2D();
        boundaryElement->SetBoundaryConditionType(boundaryConditionType);
    }

    myStructure.NodeBuildGlobalDofs();

    myStructure.Info();
    myStructure.CalculateMaximumIndependentSets();

    // apply some displacements
    int allNodes = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNodeCoordinateRange(allNodes, 0, -1, lX + 1);
    auto nodeIds = myStructure.GroupGetMemberIds(allNodes);
    for (int i = 0; i < nodeIds.GetNumRows(); ++i)
    {
        NuTo::NodeBase* node = myStructure.NodeGetNodePtr(nodeIds.GetValue(i));
        auto disps = node->GetCoordinates2D() / lX * BoundaryDisplacement;
        node->SetDisplacements2D(disps);

        if (node->GetNumNonlocalEqStrain() > 0)
            node->SetNonlocalEqStrain(disps.at(0, 0) / 10);
    }

    bool globalStiffnessCorrect = myStructure.CheckCoefficientMatrix_0(1.e-8, true);
    bool elementStiffnessCorrect = myStructure.ElementCheckCoefficientMatrix_0(1.e-8);

    if (not globalStiffnessCorrect)
        throw NuTo::Exception("global stiffness matrix incorrect!");
    if (not elementStiffnessCorrect)
        throw NuTo::Exception("element stiffness matrices incorrect!");

#ifdef ENABLE_VISUALIZE
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentNonlocalEqStrain();

    std::string resultDir = "./ResultsGradientDamage";
    boost::filesystem::create_directory(resultDir);
    myStructure.ExportVtkDataFileElements(resultDir+"/Elements.vtu", true);
    myStructure.ExportVtkDataFileNodes(resultDir+"/Nodes.vtu", true);

#endif

}


void GradientDamage3D()
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    double lX = 3, lY = 4, lZ = 5;
    NuTo::FullVector<int, Eigen::Dynamic> nodeIds(8);
    nodeIds[0] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({ 0, 0, 0}));
    nodeIds[1] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({lX, 0, 0}));
    nodeIds[2] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({lX,lY, 0}));
    nodeIds[3] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({ 0,lY, 0}));
    nodeIds[4] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({ 0, 0,lZ}));
    nodeIds[5] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({lX, 0,lZ}));
    nodeIds[6] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({lX,lY,lZ}));
    nodeIds[7] = myStructure.NodeCreate(NuTo::FullVector<double, 3> ({ 0,lY,lZ}));

    int myInterpolationType = myStructure.InterpolationTypeCreate("Brick3D");
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    myStructure.ElementCreate(myInterpolationType, nodeIds, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);

    myStructure.ElementTotalConvertToInterpolationType();

    myStructure.SetVerboseLevel(10);
    myStructure.SetShowTime(true);
    myStructure.ElementTotalConvertToInterpolationType();

    int mySection = myStructure.SectionCreate("VOLUME");
    int myConstitutiveLaw = SetConstitutiveLaw(myStructure);


    myStructure.ElementTotalSetConstitutiveLaw(myConstitutiveLaw);
    myStructure.ElementTotalSetSection(mySection);

    myStructure.NodeBuildGlobalDofs();

    myStructure.Info();
    myStructure.CalculateMaximumIndependentSets();


    // apply some displacements
    double boundaryDisplacement = 1.;

    int allNodes = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNodeCoordinateRange(allNodes, 0, -1, lX + 1);
    nodeIds = myStructure.GroupGetMemberIds(allNodes);
    for (int i = 0; i < nodeIds.GetNumRows(); ++i)
    {
        NuTo::NodeBase* node = myStructure.NodeGetNodePtr(nodeIds.GetValue(i));
        auto disps = node->GetCoordinates3D() / lX * boundaryDisplacement;
        node->SetDisplacements3D(disps);

        if (node->GetNumNonlocalEqStrain() > 0)
            node->SetNonlocalEqStrain(disps.at(0, 0) / 10);
    }

    bool globalStiffnessCorrect = myStructure.CheckCoefficientMatrix_0(1.e-8, true);
    bool elementStiffnessCorrect = myStructure.ElementCheckCoefficientMatrix_0(1.e-8);

    if (not globalStiffnessCorrect)
        throw NuTo::Exception("global stiffness matrix incorrect!");
    if (not elementStiffnessCorrect)
        throw NuTo::Exception("element stiffness matrices incorrect!");

#ifdef ENABLE_VISUALIZE
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentNonlocalEqStrain();

    std::string resultDir = "./ResultsGradientDamage";
    boost::filesystem::create_directory(resultDir);
    myStructure.ExportVtkDataFileElements(resultDir+"/Elements.vtu", true);
    myStructure.ExportVtkDataFileNodes(resultDir+"/Nodes.vtu", true);

#endif

}


void GroupRemoveNodesWithoutDisplacements(NuTo::Structure& rStructure, int rGroupNodeId)
{
    auto ids = rStructure.GroupGetMemberIds(rGroupNodeId);
    for (int i = 0; i < ids.GetNumRows(); ++i)
    {
        int nodeId = ids.GetValue(i);
        NuTo::NodeBase* node = rStructure.NodeGetNodePtr(nodeId);
        if (node->GetNumDisplacements() == 0)
        {
            NuTo::GroupBase* group = rStructure.GroupGetGroupPtr(rGroupNodeId);
            group->RemoveMember(nodeId);
        }
    }
}

void Check1D2D()
{
    double lx = 100;
    double ly = 2;
    double lz = 2;

    double numElements = 20;

    NuTo::Structure myStructure1D(1);
    NuTo::Structure myStructure2D(2);

    int interpolationType1D = myStructure1D.InterpolationTypeCreate(NuTo::Interpolation::TRUSS1D);
    myStructure1D.InterpolationTypeAdd(interpolationType1D, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure1D.InterpolationTypeAdd(interpolationType1D, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT3);
    myStructure1D.InterpolationTypeAdd(interpolationType1D, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);
    myStructure1D.InterpolationTypeSetIntegrationType(interpolationType1D, NuTo::IntegrationType::IntegrationType1D2NGauss3Ip, NuTo::IpData::STATICDATA);

    int interpolationType2D = myStructure2D.InterpolationTypeCreate(NuTo::Interpolation::QUAD2D);
    myStructure2D.InterpolationTypeAdd(interpolationType2D, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure2D.InterpolationTypeAdd(interpolationType2D, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure2D.InterpolationTypeAdd(interpolationType2D, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);
    myStructure2D.InterpolationTypeSetIntegrationType(interpolationType2D, NuTo::IntegrationType::IntegrationType2D4NGauss9Ip, NuTo::IpData::STATICDATA);

    // create nodes 1D
    int numNodes = numElements + 1; // nodes for nonlocal strain/coordinates
    double lengthElement = lx / numElements;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(1);
    for (int node = 0; node < numNodes; node++)
    {
        nodeCoordinates(0) = node * lengthElement; // two nodes per element
        myStructure1D.NodeCreate(node, nodeCoordinates);
    }
    NuTo::FullVector<int, Eigen::Dynamic> elementIncidence(2);
    for (int iElement = 0; iElement < numElements; iElement++)
    {
        elementIncidence(0) = iElement;
        elementIncidence(1) = iElement + 1;
        myStructure1D.ElementCreate(interpolationType1D, elementIncidence, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);
    }

    int nodeNum = 0;
    int numElementsY = 3;
    for (int countY = 0; countY < numElementsY + 1; countY++)
    {
        for (int countX = 0; countX < numNodes; countX++)
        {
            NuTo::FullVector<double, Eigen::Dynamic> coordinates(2);
            coordinates(0) = countX * lx / numElements;
            coordinates(1) = countY * ly / numElementsY;
            myStructure2D.NodeCreate(nodeNum, coordinates);
            nodeNum++;
        }
    }

    //create elements
    for (int countY = 0; countY < numElementsY; countY++)
    {
        for (int countX = 0; countX < numElements; countX++)
        {
            NuTo::FullVector<int, Eigen::Dynamic> nodes(4);
            nodes(0) = countX + countY * numNodes;
            nodes(1) = countX + 1 + countY * numNodes;
            nodes(2) = countX + 1 + (countY + 1) * numNodes;
            nodes(3) = countX + (countY + 1) * numNodes;
            myStructure2D.ElementCreate(interpolationType2D, nodes, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::STATICDATA);
        }
    }

    int myConstitutiveLaw1D = SetConstitutiveLaw(myStructure1D);
    int myConstitutiveLaw2D = SetConstitutiveLaw(myStructure2D);

    myStructure1D.ElementTotalConvertToInterpolationType();
    myStructure2D.ElementTotalConvertToInterpolationType();

    int mySection1D = myStructure1D.SectionCreate("Truss");
    int mySection2D = myStructure2D.SectionCreate("Plane_Stress");
    int mySection1Dr = myStructure1D.SectionCreate("Truss");
    int mySection2Dr = myStructure2D.SectionCreate("Plane_Stress");

    double alpha = 0.5;

    myStructure1D.SectionSetArea(mySection1D, lz * ly);
    myStructure2D.SectionSetThickness(mySection2D, lz);

    myStructure1D.SectionSetArea(mySection1Dr, lz * ly * alpha);
    myStructure2D.SectionSetThickness(mySection2Dr, lz * alpha);

    myStructure1D.ElementTotalSetSection(mySection1D);
    myStructure1D.ElementTotalSetConstitutiveLaw(myConstitutiveLaw1D);

    myStructure2D.ElementTotalSetSection(mySection2D);
    myStructure2D.ElementTotalSetConstitutiveLaw(myConstitutiveLaw2D);

    double kappa = 0.001;

    int damagedNodes1D = myStructure1D.GroupCreate("NODES");
    int damagedElems1D = myStructure1D.GroupCreate("ELEMENTS");
    myStructure1D.GroupAddNodeCoordinateRange(damagedNodes1D, 0, lx / 2 - 1.e-4, lx / 2 + 1.e-4);
    myStructure1D.GroupAddElementsFromNodes(damagedElems1D, damagedNodes1D, false);
    auto elemIds1D = myStructure1D.GroupGetMemberIds(damagedElems1D);
    for (int i = 0; i < elemIds1D.GetNumRows(); ++i)
    {
//        myStructure1D.ElementSetSection(elemIds1D.GetValue(i),mySection1Dr);
        NuTo::ElementBase* element = myStructure1D.ElementGetElementPtr(elemIds1D.GetValue(i));
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
        {
            element->GetStaticData(i)->AsGradientDamage1D()->SetKappa(kappa);
            std::cout << element->GetIntegrationType()->GetStrIdentifier() << std::endl;
            std::cout << element->GetInterpolationType()->GetCurrentIntegrationType()->GetStrIdentifier() << std::endl;
        }
    }

    int damagedNodes2D = myStructure2D.GroupCreate("NODES");
    int damagedElems2D = myStructure2D.GroupCreate("ELEMENTS");
    myStructure2D.GroupAddNodeCoordinateRange(damagedNodes2D, 0, lx / 2 - 1.e-4, lx / 2 + 1.e-4);
    myStructure2D.GroupAddElementsFromNodes(damagedElems2D, damagedNodes2D, false);
    auto elemIds2D = myStructure2D.GroupGetMemberIds(damagedElems2D);
    for (int i = 0; i < elemIds2D.GetNumRows(); ++i)
    {
//        myStructure2D.ElementSetSection(elemIds2D.GetValue(i),mySection2Dr);
        NuTo::ElementBase* element = myStructure2D.ElementGetElementPtr(elemIds2D.GetValue(i));
        for (int i = 0; i < element->GetNumIntegrationPoints(); ++i)
        {
            element->GetStaticData(i)->AsGradientDamage1D()->SetKappa(kappa);
            std::cout << element->GetIntegrationType()->GetStrIdentifier() << std::endl;
            std::cout << element->GetInterpolationType()->GetCurrentIntegrationType()->GetStrIdentifier() << std::endl;
        }
    }

//    return;

    int leftNodes1D = myStructure1D.GroupCreate("NODES");
    myStructure1D.GroupAddNodeCoordinateRange(leftNodes1D, 0, -1.e-4, 1.e-4);

    int leftNodes2D = myStructure2D.GroupCreate("NODES");
    myStructure2D.GroupAddNodeCoordinateRange(leftNodes2D, 0, -1.e-4, 1.e-4);

    int rightNodes1D = myStructure1D.GroupCreate("NODES");
    myStructure1D.GroupAddNodeCoordinateRange(rightNodes1D, 0, lx - 1.e-4, lx + 1.e-4);

    int rightNodes2D = myStructure2D.GroupCreate("NODES");
    myStructure2D.GroupAddNodeCoordinateRange(rightNodes2D, 0, lx - 1.e-4, lx + 1.e-4);

    NuTo::FullVector<double, Eigen::Dynamic> directionX(2);
    NuTo::FullVector<double, Eigen::Dynamic> directionY(2);

    directionX << 1, 0;
    directionY << 0, 1;

    GroupRemoveNodesWithoutDisplacements(myStructure1D, leftNodes1D);
    GroupRemoveNodesWithoutDisplacements(myStructure1D, rightNodes1D);
    GroupRemoveNodesWithoutDisplacements(myStructure2D, leftNodes2D);
    GroupRemoveNodesWithoutDisplacements(myStructure2D, rightNodes2D);

    myStructure1D.ConstraintLinearSetDisplacementNodeGroup(leftNodes1D, directionX, 0.0);
    myStructure2D.ConstraintLinearSetDisplacementNodeGroup(leftNodes2D, directionX, 0.0);

    int bc1D = myStructure1D.ConstraintLinearSetDisplacementNodeGroup(rightNodes1D, directionX, 0.0);
    int bc2D = myStructure2D.ConstraintLinearSetDisplacementNodeGroup(rightNodes2D, directionX, 0.0);

    myStructure2D.ConstraintLinearSetDisplacementNode(0, directionY, 0.);

#ifdef ENABLE_VISUALIZE
    myStructure1D.AddVisualizationComponentDisplacements();
    myStructure1D.AddVisualizationComponentEngineeringStrain();
    myStructure1D.AddVisualizationComponentEngineeringStress();
    myStructure1D.AddVisualizationComponentSection();
    myStructure1D.AddVisualizationComponentDamage();
    myStructure1D.AddVisualizationComponentNonlocalEqStrain();

    myStructure2D.AddVisualizationComponentDisplacements();
    myStructure2D.AddVisualizationComponentEngineeringStrain();
    myStructure2D.AddVisualizationComponentEngineeringStress();
    myStructure2D.AddVisualizationComponentSection();
    myStructure2D.AddVisualizationComponentDamage();
    myStructure2D.AddVisualizationComponentNonlocalEqStrain();

    myStructure1D.NodeBuildGlobalDofs();
    myStructure2D.NodeBuildGlobalDofs();
    myStructure1D.CalculateMaximumIndependentSets();
    myStructure2D.CalculateMaximumIndependentSets();
#endif //ENABLE_VISUALIZE

    NuTo::NewmarkDirect myIntegrationScheme1D(&myStructure1D);
    NuTo::NewmarkDirect myIntegrationScheme2D(&myStructure2D);

    double simulationTime = 1;

    double dispEnd = 0.1;

    double timeSol = 1;

    int numLoadSteps = 100 * timeSol;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> timeDepDisp(2, 2);
    timeDepDisp << 0, 0, simulationTime, dispEnd;

    myIntegrationScheme1D.SetTimeDependentConstraint(bc1D, timeDepDisp);
    myIntegrationScheme1D.SetMaxTimeStep(simulationTime / numLoadSteps);
    myIntegrationScheme1D.SetAutomaticTimeStepping(true);
    myIntegrationScheme1D.SetTimeStep(1e-1 * simulationTime / numLoadSteps);
    myIntegrationScheme1D.SetMinTimeStep(1e-6 * myIntegrationScheme1D.GetMaxTimeStep());
    myIntegrationScheme1D.SetToleranceForce(1e-8);
    myIntegrationScheme1D.SetMaxNumIterations(12);
    myIntegrationScheme1D.SetCheckCoefficientMatrix(false);
    myIntegrationScheme1D.SetVisualizeResidual(false);
    myIntegrationScheme1D.SetPerformLineSearch(false);

    myIntegrationScheme2D.SetTimeDependentConstraint(bc2D, timeDepDisp);
    myIntegrationScheme2D.SetMaxTimeStep(simulationTime / numLoadSteps);
    myIntegrationScheme2D.SetAutomaticTimeStepping(true);
    myIntegrationScheme2D.SetTimeStep(1e-1 * simulationTime / numLoadSteps);
    myIntegrationScheme2D.SetMinTimeStep(1e-6 * myIntegrationScheme1D.GetMaxTimeStep());
    myIntegrationScheme2D.SetToleranceForce(1e-6);
    myIntegrationScheme2D.SetMaxNumIterations(12);
    myIntegrationScheme2D.SetCheckCoefficientMatrix(false);
    myIntegrationScheme2D.SetVisualizeResidual(false);
    myIntegrationScheme2D.SetPerformLineSearch(false);

    bool deleteDirectory = true;
    std::string calcDir1D = "/home/ttitsche/tmpdir/Gradient1D";
    std::string calcDir2D = "/home/ttitsche/tmpdir/Gradient2D";

    myIntegrationScheme1D.SetResultDirectory(calcDir1D, deleteDirectory);
    myIntegrationScheme2D.SetResultDirectory(calcDir2D, deleteDirectory);

    myIntegrationScheme1D.Solve(simulationTime * timeSol);
//    myIntegrationScheme2D.Solve(simulationTime * timeSol);

//    myStructure1D.CheckCoefficientMatrix_0(1.e-6, true);
//    myStructure2D.CheckCoefficientMatrix_0(1.e-6, true);

//    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> K1D, K2D;
//    NuTo::FullVector<int,Eigen::Dynamic> d1,d2;
//
//
//    myStructure1D.ElementStiffness(0,K1D, d1, d2);
//    myStructure2D.ElementStiffness(0,K2D, d1, d2);
//
//
    std::cout << "NodeInfo 1D: " << std::endl;
    myStructure1D.NodeInfo(10);
    myStructure1D.ElementInfo(10);

    // get element stresses
    int allElements = myStructure1D.GroupCreate("Elements");
    myStructure1D.GroupAddElementFromType(allElements, interpolationType1D);

    NuTo::FullVector<int, Eigen::Dynamic> elementIds = myStructure1D.GroupGetMemberIds(allElements);
    for (int iElement = 0; iElement < elementIds.GetNumRows(); ++iElement)
    {
        int elementId = elementIds(iElement);
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> stress;
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> damage;
        myStructure1D.ElementGetEngineeringStress(elementId, stress);
        myStructure1D.ElementGetDamage(elementId, damage);
        for (int iIP = 0; iIP < stress.GetNumColumns(); ++iIP)
        {
            double numericStress = stress(0, iIP);
            std::cout << "numeric stress in element " << elementId << " at IP " << iIP << ": " << numericStress << std::endl;
            std::cout << "        damage in element " << elementId << " at IP " << iIP << ": " << damage.GetValue(0, iIP) << std::endl;
        }
    }

    NuTo::FullVector<double, Eigen::Dynamic> internalGradient;

    myStructure1D.BuildGlobalGradientInternalPotentialVector(internalGradient);
    internalGradient.Info(10, 3, true);

    int elementId = 7;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> stress, strain, damage;
    myStructure1D.ElementGetEngineeringStress(elementId, stress);
    myStructure1D.ElementGetDamage(elementId, damage);
    for (int iIP = 0; iIP < stress.GetNumColumns(); ++iIP)
    {
        std::cout << "strain at IP " << iIP << ": " << strain.GetValue(0, iIP) << std::endl;
        std::cout << "damage at IP " << iIP << ": " << damage.GetValue(0, iIP) << std::endl;

        std::cout << "(1-w)*E*strain" << (1 - damage.GetValue(0, iIP)) * strain.GetValue(0, iIP) * 30000 << std::endl;

        std::cout << "stress at IP " << iIP << ": " << stress.GetValue(0, iIP) << std::endl << std::endl;
    }

//
//    std::cout << "NodeInfo 2D: " << std::endl;
//    myStructure2D.NodeInfo(10);
//    myStructure2D.ElementInfo(10);
//
//
//
//    std::cout << "K1D: " << std::endl;
//    K1D.Info(10,3,true);
//
//    std::cout << "K2D: " << std::endl;
//    K2D.Info(10,3,true);

}

int main()
{

//    try
//    {
        CheckLocalEqStrainDerivatives3D();
        CheckLocalEqStrainDerivatives();
        CheckDamageLaws();
        GradientDamage1D();
        GradientDamage2D();
        GradientDamage3D();
//        Check1D2D();

//    } catch (NuTo::MechanicsException& e)
//    {
//        std::cout << e.ErrorMessage();
//        return -1;
//    } catch (...)
//    {
//        std::cout << "Something else went wrong." << std::endl;
//        return -1;
//    }

    std::cout << std::endl;
    std::cout << "#####################################" << std::endl;
    std::cout << "##  Congratulations! Everything    ##" << std::endl;
    std::cout << "##   went better than expected!    ##" << std::endl;
    std::cout << "#####################################" << std::endl;

    return 0;

}
