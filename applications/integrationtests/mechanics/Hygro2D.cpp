#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/constitutive/laws/PorousMediaAdapter.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/elements/ContinuumBoundaryElement.h"

using namespace NuTo;

int main()
{
    Structure structure(2);

    int group, interpolationType;
    std::tie(group, interpolationType) = MeshGenerator::Grid(structure, {1.0, 1.0}, {1, 1});
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::CAPILLARY_PRESSURE,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, Node::eDof::GAS_PRESSURE,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    PorousMediaAdapter porousMedium(0.0512, 18.6, 2.27, 0.1);

    for (int id : structure.GroupGetMemberIds(group))
        structure.ElementSetConstitutiveLaw(structure.ElementGetElementPtr(id), &porousMedium);

    auto someSection = SectionPlane::Create(1.0, false);
    structure.ElementTotalSetSection(someSection);

    NodeBase* nodePtr = structure.NodeGetNodePtr(0);
    nodePtr->Set(Node::eDof::CAPILLARY_PRESSURE, 0, 20.0);
    nodePtr->Set(Node::eDof::GAS_PRESSURE, 0, 1.0);
    nodePtr = structure.NodeGetNodePtr(1);
    nodePtr->Set(Node::eDof::CAPILLARY_PRESSURE, 0, 20.0);
    nodePtr->Set(Node::eDof::GAS_PRESSURE, 0, 1.0);

    auto hessian0 = structure.BuildGlobalHessian0();
    auto hessian1 = structure.BuildGlobalHessian1();
    auto f_int = structure.BuildGlobalInternalGradient();

    std::cout << hessian0 << std::endl;
    std::cout << hessian1 << std::endl;
    std::cout << f_int << std::endl;
    return 0;
}
