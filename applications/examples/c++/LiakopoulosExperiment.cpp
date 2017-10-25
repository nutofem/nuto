#include <iostream>
#include <cfenv>
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/constitutive/laws/PorousMediaAdapter.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/elements/ContinuumBoundaryElement.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/DirectionEnum.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
#include "mechanics/PDEs/Sand.h"


using namespace NuTo;

int main()
{
    feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
    Structure structure(2);
    structure.SetNumTimeDerivatives(1);

    int group, interpolationType;
    std::tie(group, interpolationType) = MeshGenerator::Grid(structure, {0.1, 1.0}, {5, 5});
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::CAPILLARY_PRESSURE,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::GAS_PRESSURE,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    Sand sand;
    std::cout << sand.GasRelativePermeability(0.0) << std::endl;
    PorousMediaAdapter porousMedium(sand);

    for (int id : structure.GroupGetMemberIds(group))
        structure.ElementSetConstitutiveLaw(structure.ElementGetElementPtr(id), &porousMedium);

    auto someSection = SectionPlane::Create(1.0, false);
    structure.ElementTotalSetSection(someSection);

    for (int i = 0; i < structure.GetNumNodes(); ++i)
    {
        NodeBase* nodePtr = structure.NodeGetNodePtr(i);
        nodePtr->Set(Node::eDof::GAS_PRESSURE, 0, 0.1);
    }

    const double p_atm = 0.1;
    auto& bottomNodes = structure.GroupGetNodeCoordinateRange(eDirection::Y, - 1e-6, 1e-6);
    structure.Constraints().Add(Node::eDof::GAS_PRESSURE, Constraint::Value(bottomNodes, p_atm));

    auto& topNodes = structure.GroupGetNodeCoordinateRange(eDirection::Y, 1.0 - 1e-6, 1.0 + 1e-6);
    structure.Constraints().Add(Node::eDof::GAS_PRESSURE, Constraint::Value(topNodes, p_atm));

    auto& leftNodes = structure.GroupGetNodeCoordinateRange(eDirection::X, -1e-6, 1e-6);
    structure.Constraints().Add(Node::eDof::CAPILLARY_PRESSURE, Constraint::Value(leftNodes));

    auto& rightNodes = structure.GroupGetNodeCoordinateRange(eDirection::X, 0.1 - 1e-6, 0.1 + 1e-6);
    structure.Constraints().Add(Node::eDof::CAPILLARY_PRESSURE, Constraint::Value(rightNodes));

    auto hessian0 = structure.BuildGlobalHessian0();
    auto hessian1 = structure.BuildGlobalHessian0();

    std::cout << hessian0 << std::endl;
    std::cout << hessian1 << std::endl;

    NewmarkDirect newmark(&structure);

    newmark.PostProcessing().SetResultDirectory("./results", true);
    newmark.SetTimeStep(60.0);
    newmark.Solve(120.0*60.0);
    return 0;
}
