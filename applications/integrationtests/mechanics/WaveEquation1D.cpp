#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

#include "mechanics/mesh/MeshGenerator.h"

#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"

#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/Group.h"

#include "mechanics/timeIntegration/RK4.h"
#include "mechanics/timeIntegration/StructureExplicit2ndOrder.h"
#include "mechanics/timeIntegration/TimeControl.h"

#include "mechanics/structures/Assembler.h"

#include <boost/filesystem.hpp>

#include "BoostUnitTest.h"

/* A dynamic linear elastic bar with prescribed displacement at the left end
 * slowly rising from 0 to 1. This disturbance will travel with speed v = sqrt(E/rho)
 * to the right end. The calculation is done with an explicit method and the governing equation
 * (wave equation) is transformed to a first order system using StructureExplicit2ndOrder.
*/
class TestStructure : public NuTo::Structure
{
public:
    double E;
    double rho;

    TestStructure()
        : NuTo::Structure(1)
    {
        // parameters
        E = 1.;
        rho = 1.;
        double lX = 1.;
        int numElm = 40;
        SetShowTime(false);
        SetVerboseLevel(0);
        SetNumTimeDerivatives(2);

        std::pair<int, int> meshInfo = NuTo::MeshGenerator::Grid(*this, {0.}, {lX}, {numElm});

        InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS,
                             NuTo::Interpolation::eTypeOrder::LOBATTO2);

        ElementTotalConvertToInterpolationType();

        int myLaw = ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        ConstitutiveLawSetParameterDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
        ConstitutiveLawSetParameterDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.);
        ConstitutiveLawSetParameterDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, rho);
        ElementTotalSetConstitutiveLaw(myLaw);
        NodeBuildGlobalDofs();

        ElementTotalSetSection(NuTo::SectionTruss::Create(1.0));

        AddVisualizationComponent(meshInfo.first, NuTo::eVisualizeWhat::DISPLACEMENTS);

        NuTo::NodeBase& nodeLeft = NodeGetAtCoordinate(0.);
        Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                          NuTo::Constraint::Value(nodeLeft, [&](double t) { return DirichletBoundaryCondition(t); }));

        NodeBuildGlobalDofs();
    }

    double DirichletBoundaryCondition(double t)
    {
        double tau = 0.2;
        if (t < 0)
            return 0.;
        if (t <= tau)
            return (0.5 - cos(M_PI * t / tau) / 2);
        return 1.;
    }

    double ExpectedResult(double x, double t)
    {
        return DirichletBoundaryCondition(t - x / sqrt(E / rho));
    }
};

void Solve(TestStructure& s, double timeStep, int numSteps)
{
    // delete result directory if it exists and create it new
    boost::filesystem::path resultDirectory =
            boost::filesystem::initial_path().string() + std::string("/WaveEquation1D_Results");

    NuTo::TimeIntegration::RK4<NuTo::TimeIntegration::StructureStateExplicit2ndOrder> ti;
    NuTo::TimeControl timeControl;
    NuTo::PostProcessor myPost(s, timeControl);
    myPost.SetResultDirectory(resultDirectory.string(), true);
    myPost.SetMinTimeStepPlot(numSteps * timeStep * 0.01);
    myPost.AddResultTime("Time");

    NuTo::StructureOutputBlockVector dummy(s.GetDofStatus());

    // Set up initial value
    NuTo::TimeIntegration::StructureStateExplicit2ndOrder x(s.NodeExtractDofValues(0).J, s.NodeExtractDofValues(1).J);

    // Generate ODESystem right hand side
    NuTo::TimeIntegration::StructureRhsExplicit2ndOrder eqSystem(s);

    for (int i = 0; i < numSteps; i++)
    {
        double t = i * timeStep;
        x = ti.DoStep(eqSystem, x, t, timeStep);
        //
        timeControl.SetCurrentTime(t);
        myPost.PostProcess(dummy);
        myPost.SetCallback([&](const NuTo::StructureBase& str, const NuTo::TimeControl tCtrl) {
            // ----------------------------------------------------
            // Loop over all nodes and compare with expected result
            // ----------------------------------------------------
            double maxerror = 0.;

            for (int curNode : s.GroupGetMemberIds(s.GroupGetNodesTotal()))
            {
                Eigen::VectorXd coordinates(1);
                Eigen::VectorXd displ(1);
                s.NodeGetCoordinates(curNode, coordinates);
                s.NodeGetDisplacements(curNode, displ);

                double expected = s.ExpectedResult(coordinates(0), tCtrl.GetCurrentTime());
                double computed = displ(0);

                maxerror = std::max(std::abs(computed - expected), maxerror);
            }
            BOOST_CHECK_SMALL(maxerror, 0.01);
        });
    }
}

BOOST_AUTO_TEST_CASE(TimeDependentDirichletBoundary)
{
    int numSteps = 300;
    double timeStep = 0.001;

    TestStructure s;
    Solve(s, timeStep, numSteps);
}
