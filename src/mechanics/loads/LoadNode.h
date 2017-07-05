#pragma once

#include "mechanics/loads/LoadBase.h"

namespace NuTo
{
class NodeBase;

//! @brief Abstract class for all constraints applied to a single node
class LoadNode : public LoadBase
{

public:
    //! @brief Constructor
    LoadNode(const NodeBase* rNode);

protected:
    LoadNode() = default;
    const NodeBase* mNode;
};
}//namespace NuTo


