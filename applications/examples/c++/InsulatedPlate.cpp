// This is a finite element solution of Example 3.4 from "Jiji: Heat Conduction, 3rd Edition".
// Consider a retangular plate of dimension L x H.
// On the left, T = 0, on the right T = T_0 * y.
// Top and bottom are insulated (no heat flux).
#include <iostream>
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/tools/MeshGenerator.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"

int main()
{
    double L = 2.0;
    double H = 1.0;
    double T_0 = 1.0;
    NuTo::Structure structure(2);
    structure.SetNumTimeDerivatives(0);

    auto section = structure.SectionCreate(NuTo::eSectionType::PLANE_STRAIN);
    structure.SectionSetThickness(section, 1.0);
    auto law = structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(
            law, NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.0);

    auto interpolationType = structure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(
            interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(
            interpolationType, NuTo::Node::eDof::TEMPERATURE, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(
            interpolationType, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);

    NuTo::MeshGenerator::MeshRectangularPlane(structure, section, law, interpolationType, {10, 10}, {L, H});

    structure.ElementTotalConvertToInterpolationType();

    auto nodesLeft = structure.GroupCreate(NuTo::eGroupId::Nodes);
    auto nodesRight = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(nodesLeft, 0, 0.0, 0.0);
    structure.GroupAddNodeCoordinateRange(nodesRight, 0, L, L);

    structure.ConstraintLinearSetTemperatureNodeGroup(nodesLeft, 0.0);

    auto nodeGroup = dynamic_cast<NuTo::Group<NuTo::NodeBase>*>(structure.GroupGetGroupPtr(nodesRight));
    for (auto& node : *nodeGroup)
    {
        auto coordinate = node.second->Get(NuTo::Node::eDof::COORDINATES);
        structure.ConstraintLinearSetTemperatureNode(node.second, T_0 * coordinate[1]);
    }

    structure.SolveGlobalSystemStaticElastic();

    auto allElements = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(allElements);
    structure.AddVisualizationComponent(allElements, NuTo::eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(allElements, NuTo::eVisualizeWhat::HEAT_FLUX);
    structure.ExportVtkDataFileElements("InsulatedPlate.vtk");
    return 0;
}
