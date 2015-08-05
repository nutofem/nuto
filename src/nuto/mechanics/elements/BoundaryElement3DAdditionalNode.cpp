#include "BoundaryElement3DAdditionalNode.h"


NuTo::BoundaryElement3DAdditionalNode::BoundaryElement3DAdditionalNode(const ElementBase* rBaseElement, int rSurfaceEdge, NodeBase *rNodeDependent)
    : BoundaryElement3D::BoundaryElement3D(rBaseElement, rSurfaceEdge),
      mNodeDependent(rNodeDependent)
{

}

//! @brief Gets the additional node of an boundary element, if it has one
//! @return Additional boundary node
NuTo::NodeBase *NuTo::BoundaryElement3DAdditionalNode::GetAdditionalBoundaryNode() const
{
    return mNodeDependent;
}


