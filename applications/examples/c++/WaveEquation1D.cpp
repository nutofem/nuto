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

void Solve(NuTo::Structure& s, double timeStep, int numSteps)
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
    }
}

class TestStructure : public NuTo::Structure
{
public:
    TestStructure()
        : NuTo::Structure(1)
    {
        // parameters
        double lX = 1.;
        int numElm = 40;
        double tau = 0.2;
        SetShowTime(false);
        SetVerboseLevel(0);
        SetNumTimeDerivatives(2);

        std::pair<int, int> meshInfo = NuTo::MeshGenerator::Grid(*this, {0.}, {lX}, {numElm});

        InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS,
                             NuTo::Interpolation::eTypeOrder::LOBATTO4);

        ElementTotalConvertToInterpolationType();

        int myLaw = ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        ConstitutiveLawSetParameterDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.);
        ConstitutiveLawSetParameterDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.);
        ConstitutiveLawSetParameterDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, 1.);
        ElementTotalSetConstitutiveLaw(myLaw);
        NodeBuildGlobalDofs();

        ElementTotalSetSection(NuTo::SectionTruss::Create(1.0));

        AddVisualizationComponent(meshInfo.first, NuTo::eVisualizeWhat::DISPLACEMENTS);

        NuTo::NodeBase& nodeLeft = NodeGetAtCoordinate(0.);
        Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Value(nodeLeft, [tau](double t) {
                              if (t <= tau)
                                  return (0.5 - cos(M_PI * t / tau) / 2);
                              return 1.;
                          }));

        std::cout << "Solve " << std::endl;
        NodeBuildGlobalDofs();
    }
};

int main(int argc, char* argv[])
{
    int numSteps = 1000;
    double timeStep = 0.001;

    TestStructure s;
    Solve(s, timeStep, numSteps);
}
