#include "BoostUnitTest.h"
#include <boost/filesystem.hpp>

#include "math/EigenCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/VelocityVerlet.h"
#include "mechanics/timeIntegration/RungeKutta4.h"
#include "mechanics/timeIntegration/RungeKutta2.h"
#include "mechanics/timeIntegration/RungeKutta3.h"
#include "mechanics/timeIntegration/RungeKutta38.h"
#include "mechanics/timeIntegration/NystroemQinZhu.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/constraints/ConstraintCompanion.h"

// Test of explicit time integration schemes (Kunge Kutta 4, Velocity Verlet and Nystroem Qin Zhu)
// Dynamic simulation of a 1D elastic bar under shock load
void Run(NuTo::Structure& s, NuTo::TimeIntegrationBase& timeIntegrationScheme)
{
    std::string resultDir = boost::filesystem::initial_path().string() + std::string("/ResultsExplicitTimeIntegration");
    // delete result directory
    if (boost::filesystem::exists(resultDir)) // does p actually exist?
    {
        if (boost::filesystem::is_directory(resultDir)) // is p a directory?
        {
            boost::filesystem::remove_all(resultDir);
        }
    }
    // create result directory
    boost::filesystem::create_directory(resultDir);

    // INPUT
    // geometry parameter
    double mL = 20; // length [m]
    double mA = 25e-4; // cross section [m²]

    // material parameter
    double youngsModulus = 70600e6; // Youngs modulus [N/m2] aluminium
    double poissonRatio = 0.; // Poisson ratio
    double density = 2824; // density [kg/m3]

    double THit = 0.001; // time of shock
    double mFreq = 1; // in this load case only because of time step computation
    double mAmplS = -1e6;

    // mesh parameter
    int numElements = 50; // number of elements
    double timeStep = 1e-6;

    // computation simulation time
    double waveSpeedL =
            sqrt(youngsModulus * (1. - poissonRatio) / (density * (1. - 2. * poissonRatio) * (1. + poissonRatio)));
    double simulationTime(1. * mL / waveSpeedL);
    double timePeriod(1 / mFreq);

    // set number of time derivatives to 2 (nodes have disp, vel and accelerations)
    s.SetNumTimeDerivatives(2);

    auto meshInfo = NuTo::MeshGenerator::Grid(s, {mL}, {numElements});
    s.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.ElementTotalConvertToInterpolationType();

    // create section
    auto mySection = NuTo::SectionTruss::Create(mA);
    s.ElementTotalSetSection(mySection);

    // create constitutive law
    using namespace NuTo::Constitutive;
    int myMaterial = s.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    s.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::POISSONS_RATIO, poissonRatio);
    s.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::DENSITY, density);
    s.ElementTotalSetConstitutiveLaw(myMaterial);

    // create NodeGroups
    auto& nodeLeft = s.NodeGetAtCoordinate(Eigen::Matrix<double, 1, 1>::Constant(0));
    auto& nodeRight = s.NodeGetAtCoordinate(Eigen::Matrix<double, 1, 1>::Constant(mL));

    s.CalculateMaximumIndependentSets();

    // Info
    s.SetShowTime(false);
    s.Info();

    // time step (explicit method)
    double minTimeStepAccuracy(timePeriod / 50); // to cover the external wave
    // calculate critical time step
    double criticalTimeStep(minTimeStepAccuracy);
    if (timeIntegrationScheme.HasCriticalTimeStep())
    {
        double criticalTimeStepMethod = timeIntegrationScheme.CalculateCriticalTimeStep();
        if (criticalTimeStepMethod < minTimeStepAccuracy)
        {
            criticalTimeStep = criticalTimeStepMethod;
        }
    }
    if (timeStep == 0 || timeStep > criticalTimeStep)
    {
        timeStep = criticalTimeStep;
    }
    timeIntegrationScheme.SetTimeStep(timeStep);

    // set load
    Eigen::Matrix<double, 4, 2> forceRHS;
    forceRHS(0, 0) = 0;
    forceRHS(0, 1) = mAmplS;
    forceRHS(1, 0) = THit;
    forceRHS(1, 1) = mAmplS;
    forceRHS(2, 0) = THit + timeStep;
    forceRHS(2, 1) = 0.0;
    forceRHS(3, 0) = simulationTime + 1;
    forceRHS(3, 1) = 0.0;

    // apply load
    double unitLoad = 1.;
    s.LoadCreateNodeForce(s.NodeGetId(&nodeRight), Eigen::Matrix<double, 1, 1>::Constant(1), unitLoad);
    timeIntegrationScheme.SetTimeDependentLoadCase(0, forceRHS);

    // fixed displacement at origin
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Value(nodeLeft));

    // set output
    s.SetShowTime(false);
    s.SetNumProcessors(1);
    s.SetVerboseLevel(0);

    timeIntegrationScheme.PostProcessing().AddResultTime("Time");
    timeIntegrationScheme.PostProcessing().AddResultNodeDisplacements("DisplacementsNodeRight", s.NodeGetId(&nodeRight));
    int plotElement = s.GetNumElements() / 2;
    timeIntegrationScheme.PostProcessing().AddResultElementIpData("StressCenterElement", plotElement,
                                                 NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS);

    // only plot at every 5%
    timeIntegrationScheme.PostProcessing().SetMinTimeStepPlot(simulationTime * 0.05);

    // set result directory
    bool deleteDirectory(false);
    timeIntegrationScheme.PostProcessing().SetResultDirectory(resultDir, deleteDirectory);

    s.NodeBuildGlobalDofs();

    // solve (perform explicit time integration)
    timeIntegrationScheme.SetShowTime(false);
    timeIntegrationScheme.Solve(simulationTime);

    // analytical solution
    double uanal = mAmplS * waveSpeedL * THit / (youngsModulus * mA);

    // read in result file
    boost::filesystem::path resultFile = resultDir;
    resultFile /= std::string("DisplacementsNodeRight.dat");

    Eigen::MatrixXd result = NuTo::EigenCompanion::ReadFromFile(resultFile.string());
    int EndR = result.rows();
    double difference = result(EndR - 1) - uanal;

    BOOST_CHECK_SMALL(difference, 1.e-4);
}


BOOST_AUTO_TEST_CASE(ExplicitTIRK4)
{
    NuTo::Structure s(1);
    NuTo::RungeKutta4 ti(&s);
    Run(s, ti);
}

BOOST_AUTO_TEST_CASE(ExplicitTIRK2)
{
    NuTo::Structure s(1);
    NuTo::RungeKutta2 ti(&s);
    Run(s, ti);
}

BOOST_AUTO_TEST_CASE(ExplicitTIRK3)
{
    NuTo::Structure s(1);
    NuTo::RungeKutta3 ti(&s);
    Run(s, ti);
}

BOOST_AUTO_TEST_CASE(ExplicitTIRK38)
{
    NuTo::Structure s(1);
    NuTo::RungeKutta38 ti(&s);
    Run(s, ti);
}

BOOST_AUTO_TEST_CASE(ExplcitVelocityVerlet)
{
    NuTo::Structure s(1);
    NuTo::VelocityVerlet ti(&s);
    Run(s, ti);
}

BOOST_AUTO_TEST_CASE(ExplicitNystroemQinZhu)
{
    NuTo::Structure s(1);
    NuTo::NystroemQinZhu ti(&s);
    Run(s, ti);
}
