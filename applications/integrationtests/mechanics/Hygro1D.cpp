#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/constitutive/laws/PorousMediaAdapter.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/nodes/NodeBase.h"

using namespace NuTo;

int main()
{
    Structure structure(1);

    int group, interpolationType;
    std::tie(group, interpolationType) = MeshGenerator::Grid(structure, {1.0}, {1});
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::CAPILLARY_PRESSURE,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::GAS_PRESSURE,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    PorousMediaAdapter porousMedium(0.5, 20.0, 2.0, 0.1);

    for (int id : structure.GroupGetMemberIds(group))
        structure.ElementSetConstitutiveLaw(structure.ElementGetElementPtr(id), &porousMedium);

    auto someSection = SectionTruss::Create(1.0);
    structure.ElementTotalSetSection(someSection);

    NodeBase* nodePtr = structure.NodeGetNodePtr(0);
    nodePtr->Set(Node::eDof::CAPILLARY_PRESSURE, 0, 1.0);
    nodePtr->Set(Node::eDof::GAS_PRESSURE, 0, 1.0);
    nodePtr = structure.NodeGetNodePtr(1);
    nodePtr->Set(Node::eDof::CAPILLARY_PRESSURE, 0, 1.0);
    nodePtr->Set(Node::eDof::GAS_PRESSURE, 0, 1.0);

    auto hessian0 = structure.BuildGlobalHessian0();
    auto hessian1 = structure.BuildGlobalHessian1();

    std::cout << hessian0 << std::endl;
    std::cout << hessian1 << std::endl;
    return 0;
}
