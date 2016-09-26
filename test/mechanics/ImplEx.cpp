#include <boost/filesystem.hpp>

#include "nuto/base/Timer.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/sections/SectionTruss.h"

#include "nuto/mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/timeIntegration/ImplEx.h"

#include "nuto/mechanics/elements/ContinuumBoundaryElement.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include "nuto/mechanics/constitutive/staticData/Leaf.h"


#include "nuto/mechanics/elements/IpDataStaticData.h"
#include "nuto/visualize/VisualizeEnum.h"

namespace NuToTest {
namespace ImplEx {


int SetConstitutiveLaw(NuTo::Structure& rStructure)
{
    // create a damage law
    int lawId = rStructure.ConstitutiveLawCreate("Gradient_Damage_Engineering_Stress");
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::DENSITY, 1.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 1);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 0.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    rStructure.ConstitutiveLawSetParameterDouble(lawId,NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    rStructure.ConstitutiveLawSetDamageLaw(lawId, NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

//    int myNumberConstitutiveLaw = rStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
//    rStructure.ConstitutiveLawSetDensity(myNumberConstitutiveLaw, 1.0);
//    rStructure.ConstitutiveLawSetYoungsModulus(myNumberConstitutiveLaw, 30000);
//    rStructure.ConstitutiveLawSetPoissonsRatio(myNumberConstitutiveLaw, 0.0);

    return lawId;
}

void Extrapolate()
{
    NuTo::ConstitutiveTimeStep deltaT(2);
    deltaT[0] = 1.5; // current time step
    deltaT[1] = 2.5; // previous time step

    {
        // extrapolate a double
        double X1 = 2;
        double X0 = 1;
        double result = NuTo::ConstitutiveCalculateStaticData::EulerForward(X1, X0, deltaT);
        if (std::abs(result - 2.6) > 1.e-10) throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "EulerForward<double> incorrect");
    }

    {
        // extrapolate a constitutive scalar
        NuTo::ConstitutiveScalar X1; X1[0] = 2;
        NuTo::ConstitutiveScalar X0; X0[0] = 1;
        auto result = NuTo::ConstitutiveCalculateStaticData::EulerForward(X1, X0, deltaT);
        if (std::abs(result[0] - 2.6) > 1.e-10) throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "EulerForward<ConstitutiveScalar> incorrect");
    }

    {
        // extrapolate a constitutive vector
        NuTo::ConstitutiveVector<1> X1; X1[0] = 2;
        NuTo::ConstitutiveVector<1> X0; X0[0] = 1;
        auto result = NuTo::ConstitutiveCalculateStaticData::EulerForward(X1, X0, deltaT);
        if (std::abs(result[0] - 2.6) > 1.e-10) throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "EulerForward<ConstitutiveVector> incorrect");
    }

    {
        // extrapolate a constitutive matrix
        NuTo::ConstitutiveMatrix<1,1> X1; X1(0,0) = 2;
        NuTo::ConstitutiveMatrix<1,1> X0; X0(0,0) = 1;
        auto result = NuTo::ConstitutiveCalculateStaticData::EulerForward(X1, X0, deltaT);
        if (std::abs(result(0,0) - 2.6) > 1.e-10) throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "EulerForward<ConstitutiveMatrix> incorrect");
    }




    NuTo::ConstitutiveTimeStep t(4);
    t[0] = 42;
    t.SetCurrentTimeStep(1337); // t should be [1337, 42, 0 0]
    if (not (t[0] == 1337 and t[1] == 42 and t[2] == 0 and t[3] == 0))
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "SetCurrentTimeStep incorrect");

    t.SetCurrentTimeStep(6174); // t should be [6174, 1337, 42, 0 0]
    if (not (t[0] == 6174 and t[1] == 1337 and t[2] == 42 and t[3] == 0))
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "SetCurrentTimeStep incorrect");

}


void ImplEx()
{
    NuTo::Timer timer(__FUNCTION__);
    const int numElements = 50;
    const double length = 50.;
    const double area = 10;

    NuTo::Structure s(1);
    s.SetShowTime(false);
    s.SetVerboseLevel(0);

    // create nodes
    int numNodes = numElements + 1;
    double lengthElement = length / numElements;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(1);
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        nodeCoordinates(0) = iNode * lengthElement;
        s.NodeCreate(iNode, nodeCoordinates);
    }

    int interpolationType = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSS1D);
    s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
//    s.InterpolationTypeSetIntegrationType(interpolationType, NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::STATICDATA);

    // create elements
    NuTo::FullVector<int, Eigen::Dynamic> nodes(2);
    for (int iElement = 0; iElement < numElements; iElement++)
    {
        nodes(0) = iElement;
        nodes(1) = iElement + 1;
        s.ElementCreate(interpolationType, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    }

    // create sections
    int mySection = s.SectionCreate("Truss");
    s.SectionSetArea(mySection, area);

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s));
    s.ElementTotalSetSection(mySection);

    // Set a pre-damaged spot in the middle element
    double kappa_i = 4./30000.;
    double kappa_star =  kappa_i * 3;
    auto eWeak = s.ElementGetElementPtr(numElements/2.);
    for (int i = 0; i < eWeak->GetNumIntegrationPoints(); ++i)
        dynamic_cast<NuTo::Constitutive::StaticData::Leaf<double>*>(eWeak->GetConstitutiveStaticData(i))->SetData(kappa_star);

//
//    int mySection2 = s.SectionCreate("Truss");
//    s.SectionSetArea(mySection2, area*.9);
//    s.ElementSetSection(24, mySection2);
//    s.ElementSetSection(25, mySection2);

    s.ElementGroupAllocateAdditionalStaticData(s.GroupGetElementsTotal(), 2);

    s.ConstraintLinearSetDisplacementNode(0, NuTo::FullVector<double, 1>::UnitX(), 0.0);

    int gNodeBC = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeCoordinateRange(gNodeBC, 0, length, length);
    int iNodeBC = s.GroupGetMemberIds(gNodeBC)[0];

    int bc = s.ConstraintLinearSetDisplacementNode(iNodeBC, NuTo::FullVector<double, 1>::UnitX(), 0.0);

    int visualizationGroup = s.GroupGetElementsTotal();

    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DAMAGE);

    s.CalculateMaximumIndependentSets();

    NuTo::ImplEx myIntegrationScheme(&s);

    double simulationTime = 1;
    double dispEnd = 0.1;
    int numLoadSteps = 10;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> timeDepDisp(2, 2);
    timeDepDisp << 0, 0, simulationTime, dispEnd;

    myIntegrationScheme.AddTimeDependentConstraint(bc, timeDepDisp);
    myIntegrationScheme.SetTimeStep(simulationTime / numLoadSteps);
    myIntegrationScheme.SetAutomaticTimeStepping(true);

    myIntegrationScheme.AddCalculationStep({NuTo::Node::eDof::DISPLACEMENTS});
    myIntegrationScheme.AddCalculationStep({NuTo::Node::eDof::NONLOCALEQSTRAIN});

    myIntegrationScheme.AddDofWithConstantHessian0(NuTo::Node::eDof::NONLOCALEQSTRAIN);

    myIntegrationScheme.AddResultGroupNodeForce("Force", gNodeBC);
    myIntegrationScheme.AddResultNodeDisplacements("Displ", iNodeBC);

    myIntegrationScheme.SetExtrapolationErrorThreshold(kappa_i * 0.02);

    std::string resultDir = "./ResultsImplex";
    boost::filesystem::create_directory(resultDir);
    myIntegrationScheme.SetResultDirectory(resultDir, true);

    myIntegrationScheme.Solve(simulationTime / 5.);

}

}  // namespace NuToTest
}  // namespace ImplEx

int main()
{


    try
    {
        NuTo::Timer Timer("ImplEx", true);

        NuToTest::ImplEx::Extrapolate();

        NuToTest::ImplEx::ImplEx();

    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cout << "Something else went wrong." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << std::endl;
    std::cout << "#####################################" << std::endl;
    std::cout << "##  Congratulations! Everything    ##" << std::endl;
    std::cout << "##   went better than expected!    ##" << std::endl;
    std::cout << "#####################################" << std::endl;

    return EXIT_SUCCESS;

}
