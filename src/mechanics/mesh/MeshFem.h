#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "base/Group.h"
#include "base/ValueVector.h"
#include "mechanics/DirectionEnum.h"
#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"

#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

#include <sstream>

namespace NuTo
{
class MeshFem
{
public:
    InterpolationSimple& CreateInterpolation(const InterpolationSimple& interpolation);


    //! @brief selects a coordinate at given `coords`
    //! @param coords global coordinates
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    NodeSimple& NodeAtCoordinate(Eigen::VectorXd coords, double tol = 1.e-6);

    //! @brief selects a node of type `dofType` at given `coords`
    //! @param coords global coordinates
    //! @param dofType dof type
    //! @param tol selection tolerance
    //! @return reference to the selected node, throws, if no node is found
    NodeSimple& NodeAtCoordinate(Eigen::VectorXd coords, DofType dofType, double tol = 1.e-10);

    //! @brief selects all nodes of type `dofType` where the `coord` in `direction` is within `tol`
    //! @param direction ::X, ::Y, or ::Z
    //! @param dofType dof type
    //! @param axisOffset distance of the node to the axis
    //! @param tol selection tolerance
    //! @return group with selected nodes, the group may be empty if no nodes were found
    Groups::Group<NodeSimple> NodesAtAxis(eDirection direction, DofType dofType, double axisOffset = 0.,
                                          double tol = 1.e-6);

    //! @brief selects all nodes of `dofType`
    //! @return group containing all selected nodes
    Groups::Group<NodeSimple> NodesTotal(DofType dofType);

    //! @brief selects all coordinate nodes
    //! @return group containing all selected nodes
    Groups::Group<NodeSimple> NodesTotal();

public:
    ValueVector<NodeSimple> Nodes;
    ValueVector<ElementCollectionFem> Elements;

private:
    boost::ptr_vector<InterpolationSimple> mInterpolations;
};

namespace UnitMeshFem
{

//! @brief creates a triangular mesh from (0,0) -- (1,1) with numX and numY divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @return created mesh
MeshFem CreateTriangles(int numX, int numY);

//! @brief creates a quad mesh from (0,0) -- (1,1) with numX and numY divisions
//! @param numX number of divisions in x direction
//! @param numY number of divisions in y direction
//! @return created mesh
MeshFem CreateQuads(int numX, int numY);

//! @brief transforms a mesh with a given transformation function f
//! @param rMesh mesh that is modified (return argument! pointer syntax to make it clear)
//! @param f transformation function
void Transform(MeshFem* rMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f);

} /* UnitMeshFem */
} /* NuTo */
