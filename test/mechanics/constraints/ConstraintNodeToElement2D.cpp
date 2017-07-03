//
// Created by phuschke on 2/21/17.
//

#include "BoostUnitTest.h"

#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/laws/HeatConduction.h"
#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;

constexpr int dim = 2;

int CreateNode(const Eigen::VectorXd& nodeCoordinates, NuTo::Structure &structure);
void BuildStructure(NuTo::Structure& s);

/// \brief Checks if a node can be constrained to an element
/// Creates a reference triangle and arbitrary nodes and checks whether the node coordinates are inside the element.
BOOST_AUTO_TEST_CASE(check_constraint_node_to_element)
{
    NuTo::Structure s(dim);
    BuildStructure(s);
    const int groupElementId = s.GroupGetElementsTotal();

    Eigen::VectorXd nodeCoordinates(dim);

    // this node is inside the reference tetrahedron
    nodeCoordinates = Eigen::VectorXd::Constant(dim, 0.1);
    s.ConstraintLinearEquationNodeToElementCreate(CreateNode(nodeCoordinates, s), groupElementId, eDof::DISPLACEMENTS);

    // this node is inside the reference tetrahedron
    nodeCoordinates = Eigen::VectorXd::Constant(dim, 0.);
    s.ConstraintLinearEquationNodeToElementCreate(CreateNode(nodeCoordinates, s), groupElementId, eDof::DISPLACEMENTS);

    // this node is NOT inside the reference tetrahedron
    nodeCoordinates = Eigen::VectorXd::Constant(dim, 1.);
    BOOST_REQUIRE_THROW(s.ConstraintLinearEquationNodeToElementCreate(CreateNode(nodeCoordinates, s), groupElementId, eDof::DISPLACEMENTS), NuTo::MechanicsException);

    // this node is NOT inside the reference tetrahedron but the offset lets us constrain the node anyway :-)
    Eigen::Vector3d nodeCoordinateOffset = Eigen::Vector3d::Constant(-1.);
    s.ConstraintLinearEquationNodeToElementCreate(CreateNode(nodeCoordinates, s), groupElementId, eDof::DISPLACEMENTS, 1.e-6, nodeCoordinateOffset);

    // this node is NOT inside the reference tetrahedron but the default tolerance (1.e-6) accepts it
    nodeCoordinates = Eigen::VectorXd::Constant(dim, -1.e-8);
    s.ConstraintLinearEquationNodeToElementCreate(CreateNode(nodeCoordinates, s), groupElementId, eDof::DISPLACEMENTS);

    // this node is NOT inside the reference tetrahedron and the tolerance is too small
    nodeCoordinates = Eigen::VectorXd::Constant(dim, -1.e-8);
    BOOST_REQUIRE_THROW(s.ConstraintLinearEquationNodeToElementCreate(CreateNode(nodeCoordinates, s), groupElementId, eDof::DISPLACEMENTS, 1.e-9), NuTo::MechanicsException);
}

int CreateNode(const Eigen::VectorXd& nodeCoordinates, NuTo::Structure &structure)
{
    std::set<eDof> dofSet;
    dofSet.emplace(eDof::COORDINATES);
    dofSet.emplace(eDof::DISPLACEMENTS);
    return structure.NodeCreate(nodeCoordinates, dofSet);
}

void BuildStructure(NuTo::Structure& s)
{

    s.SetShowTime(false);

    const int interpolationTypeId = s.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    s.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);

    Eigen::MatrixXd nodesCoordinates(dim,3);
    nodesCoordinates << 0, 1, 0,
                        0, 0, 1;

    auto nodeIds = s.NodesCreate(nodesCoordinates);
    s.ElementCreate(interpolationTypeId,nodeIds);
    s.ElementTotalConvertToInterpolationType();


}
