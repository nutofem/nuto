#include "mechanics/elements/ContinuumBoundaryElementConstrainedControlNode.h"

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/elements/ElementEnum.h"



template <int TDim>
NuTo::ContinuumBoundaryElementConstrainedControlNode<TDim>::ContinuumBoundaryElementConstrainedControlNode(const ContinuumElement<TDim>& rBaseElement,
                                                                                                           int rSurfaceId,
                                                                                                           NodeBase* rControlNode)
    : ContinuumBoundaryElement<TDim>::ContinuumBoundaryElement(rBaseElement,rSurfaceId),
      mControlNode(rControlNode)
{}



//! @brief returns the enum (type of the element)
//! @return enum
template <int TDim>
NuTo::Element::eElementType NuTo::ContinuumBoundaryElementConstrainedControlNode<TDim>::GetEnumType() const
{
    return Element::eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE;
}




//! @brief Gets the control node of an boundary element, if it has one
//! @return boundary control node
template <int TDim>
NuTo::NodeBase *NuTo::ContinuumBoundaryElementConstrainedControlNode<TDim>::GetBoundaryControlNode() const
{
    return mControlNode;
}



template class NuTo::ContinuumBoundaryElementConstrainedControlNode<1>;
template class NuTo::ContinuumBoundaryElementConstrainedControlNode<2>;
template class NuTo::ContinuumBoundaryElementConstrainedControlNode<3>;
