#pragma once


#include "mechanics/elements/ContinuumBoundaryElement.h"

namespace NuTo
{

template <int TDim>
class ContinuumBoundaryElementConstrainedControlNode : public ContinuumBoundaryElement<TDim>
{

public:
    ContinuumBoundaryElementConstrainedControlNode(const ContinuumElement<TDim>& rBaseElement, int rSurfaceId, NodeBase* rControlNode);

    //! @brief Gets the control node of an boundary element, if it has one
    //! @return boundary control node
    virtual NodeBase* GetBoundaryControlNode() const override;

protected:
    NodeBase* mControlNode;
};

}//namespace NuTo
