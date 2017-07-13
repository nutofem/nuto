#include "mechanics/elements/ContinuumBoundaryElementConstrainedControlNode.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/elements/ElementEnum.h"


template <int TDim>
NuTo::ContinuumBoundaryElementConstrainedControlNode<TDim>::ContinuumBoundaryElementConstrainedControlNode(
        const ContinuumElement<TDim>& rBaseElement, const IntegrationTypeBase& integrationType, int rSurfaceId,
        NodeBase* rControlNode)
    : ContinuumBoundaryElement<TDim>::ContinuumBoundaryElement(rBaseElement, integrationType, rSurfaceId)
    , mControlNode(rControlNode)
{
}


template <int TDim>
NuTo::NodeBase* NuTo::ContinuumBoundaryElementConstrainedControlNode<TDim>::GetBoundaryControlNode() const
{
    return mControlNode;
}


template class NuTo::ContinuumBoundaryElementConstrainedControlNode<1>;
template class NuTo::ContinuumBoundaryElementConstrainedControlNode<2>;
template class NuTo::ContinuumBoundaryElementConstrainedControlNode<3>;
