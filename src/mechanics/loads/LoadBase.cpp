// $Id$

#include "mechanics/loads/LoadBase.h"

//! @brief constructor
NuTo::LoadBase::LoadBase(int rLoadCase)
{
    mLoadCase = rLoadCase;
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::LoadBase)
#endif
