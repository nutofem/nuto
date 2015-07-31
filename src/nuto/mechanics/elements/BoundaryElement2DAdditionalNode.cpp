#include "BoundaryElement2DAdditionalNode.h"


NuTo::BoundaryElement2DAdditionalNode::BoundaryElement2DAdditionalNode(const ElementBase* rBaseElement, int rSurfaceEdge, NodeBase *rNodeDependent)
    : BoundaryElement2D::BoundaryElement2D(rBaseElement, rSurfaceEdge),
      mNodeDependent(rNodeDependent)
{

}

//! @brief Gets the additional node of an boundary element, if it has one
//! @return Additional boundary node
NuTo::NodeBase *NuTo::BoundaryElement2DAdditionalNode::GetAdditionalBoundaryNode() const
{
    return mNodeDependent;
}


