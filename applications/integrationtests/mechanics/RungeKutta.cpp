#include "BoostUnitTest.h"
#include <sstream>
#include "boost/filesystem.hpp"

#include "mechanics/MechanicsEnums.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/timeIntegration/RungeKutta4.h"
#include "mechanics/timeIntegration/RungeKutta2.h"
#include "mechanics/timeIntegration/RungeKutta3.h"
#include "mechanics/timeIntegration/RungeKutta38.h"
#include "mechanics/timeIntegration/RungeKuttaCashKarp.h"
#include "mechanics/timeIntegration/RungeKuttaDormandPrince.h"

#include "mechanics/mesh/MeshGenerator.h"

#include "mechanics/sections/SectionPlane.h"
/*
 *  TT:
 *  This is just a test to see if the time integration scheme compiles/runs. I highly doubt the
 *  correctness of the test setup and its results.
 */

void Run(NuTo::Structure& s, NuTo::RungeKuttaBase& rTimeIntegrationScheme)
{
    std::string resultDir = boost::filesystem::initial_path().string() + std::string("/ResultsRungeKutta");
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

    double mL = 100;
    double mH = 10;
    int numElementsX = 5;
    int numElementsY = 1;
    double density = 1; // N and mm^3


    // set number of time derivatives to 2 (nodes have disp, vel and accelerations)
    s.SetNumTimeDerivatives(2);

    // section
    double thickness = 1;
    auto mySection = NuTo::SectionPlane::Create(thickness, false);
    s.ElementTotalSetSection(mySection);

    // constitutive
    int myMatLattice = s.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    using NuTo::Constitutive::eConstitutiveParameter;
    s.ConstitutiveLawSetParameterDouble(myMatLattice, eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    s.ConstitutiveLawSetParameterDouble(myMatLattice, eConstitutiveParameter::POISSONS_RATIO, 0.0);
    s.ConstitutiveLawSetParameterDouble(myMatLattice, eConstitutiveParameter::DENSITY, density);
    s.ElementTotalSetConstitutiveLaw(myMatLattice);


    auto meshInfo = NuTo::MeshGenerator::Grid(s, {mL, mH}, {numElementsX, numElementsY},
                                              NuTo::Interpolation::eShapeType::QUAD2D);

    s.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    s.ElementGroupSetSection(meshInfo.first, mySection);
    s.ElementGroupSetConstitutiveLaw(meshInfo.first, myMatLattice);


    s.ElementTotalConvertToInterpolationType();

    // create node groups bottom boundary
    const auto& nOrigin = s.NodeGetAtCoordinate(Eigen::Vector2d::Zero());
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(nOrigin, {NuTo::eDirection::Y}));

    const auto& groupLeft = s.GroupGetNodesAtCoordinate(NuTo::eDirection::X, 0);
    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(groupLeft, {NuTo::eDirection::X}));

    int grpNodes_Right = s.GroupCreate("Nodes");
    s.GroupAddNodeCoordinateRange(grpNodes_Right, 0, mL, mL);

    // calculate maximum independent sets
    s.CalculateMaximumIndependentSets();

    // set constraints
    s.LoadCreateNodeGroupForce(0, grpNodes_Right, Eigen::Vector2d::UnitX(), 1000);

    s.Info();
    NuTo::RungeKutta4 myIntegrationScheme(&s);

    double simulationTime(1);

    // set a sinusoidal quarter wave
    double period = 5;
    Eigen::MatrixXd forceRHS(101, 2);
    for (int count = 0; count < forceRHS.rows() - 1; count++)
    {
        double t = ((double)count) / ((double)forceRHS.rows() - 2.) * 0.25 * period;
        forceRHS(count, 0) = t;
        forceRHS(count, 1) = 0.1 * sin(t / period * 2. * M_PI);
        // loadRHSFactor(count,0) = 0;
    }
    forceRHS(forceRHS.rows() - 1, 0) = forceRHS(forceRHS.rows() - 2, 0) + 1;
    forceRHS(forceRHS.rows() - 1, 1) = forceRHS(forceRHS.rows() - 2, 1);

    //        forceRHS.Info();

    rTimeIntegrationScheme.SetTimeDependentLoadCase(0, forceRHS);

    rTimeIntegrationScheme.SetMaxTimeStep(10);
    rTimeIntegrationScheme.SetMinTimeStep(0.0001 * myIntegrationScheme.GetMaxTimeStep());

    // set output during the simulation to false
    s.SetShowTime(false);

    rTimeIntegrationScheme.AddResultTime("Time");
    rTimeIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Right", grpNodes_Right);
    rTimeIntegrationScheme.AddResultElementIpData(
            "StressElement", 1, NuTo::IpData::eIpStaticDataType::ENGINEERING_STRESS); // only for test issues

    // set result directory
    bool deleteResultDirectoryFirst(true);
    rTimeIntegrationScheme.SetResultDirectory(resultDir, deleteResultDirectoryFirst);

    // solve (perform Newton raphson iteration
    rTimeIntegrationScheme.Solve(simulationTime);

    Eigen::Vector2d displacementResult =
            s.NodeGetAtCoordinate(Eigen::Vector2d(mL, mH)).Get(NuTo::Node::eDof::DISPLACEMENTS);
    BOOST_CHECK_SMALL(displacementResult.y(), 1.e-10); // uniaxial in x

    // really pointless test, no analytic solution
    BOOST_CHECK_GT(displacementResult.x(), 0);
    BOOST_CHECK_LT(displacementResult.x(), 1);
}

BOOST_AUTO_TEST_CASE(RK2)
{
    NuTo::Structure mStructure(2);
    NuTo::RungeKutta2 rk(&mStructure);
    Run(mStructure, rk);
}

BOOST_AUTO_TEST_CASE(RK3)
{
    NuTo::Structure mStructure(2);
    NuTo::RungeKutta3 rk(&mStructure);
    Run(mStructure, rk);
}

BOOST_AUTO_TEST_CASE(RK38)
{
    NuTo::Structure mStructure(2);
    NuTo::RungeKutta38 rk(&mStructure);
    Run(mStructure, rk);
}

BOOST_AUTO_TEST_CASE(RK4)
{
    NuTo::Structure mStructure(2);
    NuTo::RungeKutta4 rk(&mStructure);
    Run(mStructure, rk);
}

BOOST_AUTO_TEST_CASE(RKCashKarp)
{
    NuTo::Structure mStructure(2);
    NuTo::RungeKuttaCashKarp rk(&mStructure);
    Run(mStructure, rk);
}

BOOST_AUTO_TEST_CASE(RKDormanPrince)
{
    NuTo::Structure mStructure(2);
    NuTo::RungeKuttaDormandPrince rk(&mStructure);
    Run(mStructure, rk);
}

