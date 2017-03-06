#pragma once

#include "mechanics/nodes/DofContainer.h"

namespace NuTo
{

//! @brief dof container that is also capable of performing calculations.
//! @remark We should certainly (!) discuss the name. But I wanted it to be short and catchy :)
template <typename T>
class DofVector : public DofContainer<T>
{
public:
    void operator+=(const DofVector& rOther)
    {
        for (auto i = 0; i < rOther.mData.size(); ++i)
            this->mData[i] += rOther.mData[i];
    }
};   
} /* NuTo */  
