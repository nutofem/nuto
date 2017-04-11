#include "BoostUnitTest.h"

#include <boost/filesystem.hpp>

#include "base/Timer.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/sections/SectionTruss.h"

#include "mechanics/groups/Group.h"
#include "mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/timeIntegration/ImplEx.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/MechanicsEnums.h"

#include "mechanics/mesh/MeshGenerator.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#include "visualize/VisualizeEnum.h"


int SetConstitutiveLaw(NuTo::Structure& rStructure)
{
    // create a damage law
    using namespace NuTo::Constitutive;
    int lawId = rStructure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    rStructure.ConstitutiveLawSetParameterDouble(lawId, eConstitutiveParameter::DENSITY, 1.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId, eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    rStructure.ConstitutiveLawSetParameterDouble(lawId, eConstitutiveParameter::POISSONS_RATIO, 0.0);
    rStructure.ConstitutiveLawSetParameterDouble(lawId, eConstitutiveParameter::NONLOCAL_RADIUS, 1);
    rStructure.ConstitutiveLawSetParameterDouble(lawId, eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    rStructure.ConstitutiveLawSetParameterDouble(lawId, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    rStructure.ConstitutiveLawSetParameterDouble(lawId, eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    rStructure.ConstitutiveLawSetDamageLaw(lawId, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

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
        BOOST_CHECK_CLOSE(result, 2.6, 1.e-10);
    }

    {
        // extrapolate a constitutive scalar
        NuTo::ConstitutiveScalar X1; X1[0] = 2;
        NuTo::ConstitutiveScalar X0; X0[0] = 1;
        auto result = NuTo::ConstitutiveCalculateStaticData::EulerForward(X1, X0, deltaT);
        BOOST_CHECK_CLOSE(result[0], 2.6, 1.e-10);
    }

    {
        // extrapolate a constitutive vector
        NuTo::ConstitutiveVector<1> X1; X1[0] = 2;
        NuTo::ConstitutiveVector<1> X0; X0[0] = 1;
        auto result = NuTo::ConstitutiveCalculateStaticData::EulerForward(X1, X0, deltaT);
        BOOST_CHECK_CLOSE(result[0], 2.6, 1.e-10);
    }

    {
        // extrapolate a constitutive matrix
        NuTo::ConstitutiveMatrix<1,1> X1; X1(0,0) = 2;
        NuTo::ConstitutiveMatrix<1,1> X0; X0(0,0) = 1;
        auto result = NuTo::ConstitutiveCalculateStaticData::EulerForward(X1, X0, deltaT);
        BOOST_CHECK_CLOSE(result(0,0), 2.6, 1.e-10);
    }



    NuTo::ConstitutiveTimeStep t(4);
    t[0] = 42;
    t.SetCurrentTimeStep(1337);
    BOOST_CHECK_EQUAL(t[0], 1337);
    BOOST_CHECK_EQUAL(t[1], 42);
    BOOST_CHECK_EQUAL(t[2], 0);
    BOOST_CHECK_EQUAL(t[3], 0);

    t.SetCurrentTimeStep(6174);
    BOOST_CHECK_EQUAL(t[0], 6174);
    BOOST_CHECK_EQUAL(t[1], 1337);
    BOOST_CHECK_EQUAL(t[2], 42);
    BOOST_CHECK_EQUAL(t[3], 0);
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

    int interpolationType = NuTo::MeshGenerator::Grid(s, {length}, {numElements}).second;
    s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    // create sections
    auto section = NuTo::SectionTruss::Create(area);

    s.ElementTotalConvertToInterpolationType();
    s.ElementTotalSetConstitutiveLaw(SetConstitutiveLaw(s));
    s.ElementTotalSetSection(section);
    s.DofTypeSetIsSymmetric(NuTo::Node::eDof::DISPLACEMENTS, true);
    s.DofTypeSetIsSymmetric(NuTo::Node::eDof::NONLOCALEQSTRAIN, true);


    // Set a pre-damaged spot in the middle element
    double kappa_i = 4./30000.;
    double kappa_star =  kappa_i * 3;
    auto eWeak = s.ElementGetElementPtr(numElements/2.);
    for (int i = 0; i < eWeak->GetNumIntegrationPoints(); ++i)
        eWeak->GetIPData().GetIPConstitutiveLaw(i).GetData<NuTo::GradientDamageEngineeringStress>().SetData(kappa_star);

    s.ElementGroupAllocateAdditionalStaticData(s.GroupGetElementsTotal(), 2);

    using namespace NuTo::Constraint;
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, Component(s.NodeGetAtCoordinate(Eigen::VectorXd::Zero(1)), {NuTo::eDirection::X}));
    
    auto& nodesBC = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, length);
    double simulationTime = 1;
    double dispEnd = 0.1;
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, Component(nodesBC, {NuTo::eDirection::X}, RhsRamp(simulationTime, dispEnd)));

    int visualizationGroup = s.GroupGetElementsTotal();
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DAMAGE);

    s.CalculateMaximumIndependentSets();

    NuTo::ImplEx implex(&s);

    int numLoadSteps = 10;

    Eigen::Matrix2d timeDepDisp;
    timeDepDisp << 0, 0, simulationTime, dispEnd;

    implex.SetTimeStep(simulationTime / numLoadSteps);
    implex.SetAutomaticTimeStepping(true);

    implex.AddCalculationStep({NuTo::Node::eDof::DISPLACEMENTS});
    implex.AddCalculationStep({NuTo::Node::eDof::NONLOCALEQSTRAIN});

    implex.AddDofWithConstantHessian0(NuTo::Node::eDof::NONLOCALEQSTRAIN);

    int gNodeBC = s.GroupGetId(&nodesBC);
    int iNodeBC = s.GroupGetMemberIds(gNodeBC)[0];
    implex.AddResultGroupNodeForce("Force", gNodeBC);
    implex.AddResultNodeDisplacements("Displ", iNodeBC);

    implex.SetExtrapolationErrorThreshold(1);

    std::string resultDir = "./ResultsImplex";
    boost::filesystem::create_directory(resultDir);
    implex.SetResultDirectory(resultDir, true);

    implex.Solve(simulationTime / 5.);

}

BOOST_AUTO_TEST_CASE(ImplexExtrapolate)
{
    Extrapolate();
}

BOOST_AUTO_TEST_CASE(ImplexRun)
{
    BOOST_CHECK_NO_THROW(ImplEx());
}
