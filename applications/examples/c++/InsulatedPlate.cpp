// This is a finite element solution of Example 3.4 from "Jiji: Heat Conduction, 3rd Edition".
// Consider a retangular plate of dimension L x H.
// On the left, T = 0, on the right T = T_0 * y.
// Top and bottom are insulated (no heat flux).
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/mesh/MeshGenerator.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/visualize/VisualizeEnum.h"

using namespace NuTo;

int main()
{
    double L = 2.0;
    double H = 1.0;
    double T_0 = 1.0;
    Structure structure(2);
    structure.SetNumTimeDerivatives(0);

    auto section = SectionPlane::Create(1.0, true);

    auto law = structure.ConstitutiveLawCreate(Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(law, Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.0);

    auto meshInfo = MeshGenerator::Grid(structure, {L, H}, {10, 10});
    structure.InterpolationTypeAdd(meshInfo.second, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(meshInfo.second, eIntegrationType::IntegrationType2D4NGauss4Ip);

    structure.ElementGroupSetSection(meshInfo.first, section);
    structure.ElementGroupSetConstitutiveLaw(meshInfo.first, law);

    structure.ElementTotalConvertToInterpolationType();

    auto nodesLeft = structure.GroupGetNodeCoordinateRange(eDirection::X, 0.0, 0.0);
    auto nodesRight = structure.GroupGetNodeCoordinateRange(eDirection::X, L, L);

    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(nodesLeft));

    for (auto& node : nodesRight)
    {
        auto coordinate = node.second->Get(Node::eDof::COORDINATES);
        structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(*node.second, T_0 * coordinate[1]));
    }

    structure.SolveGlobalSystemStaticElastic();

    auto allElements = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsTotal(allElements);
    structure.AddVisualizationComponent(allElements, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(allElements, eVisualizeWhat::HEAT_FLUX);
    structure.ExportVtkDataFileElements("InsulatedPlate.vtk");
    return 0;
}
