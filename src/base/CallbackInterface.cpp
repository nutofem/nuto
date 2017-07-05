// $Id$

#include "base/CallbackInterface.h"
#include "base/Exception.h"

// serializes the class

bool NuTo::CallbackInterface::Exit(NuTo::StructureBase &)
{
    throw Exception(__PRETTY_FUNCTION__, "not implemented for this callback.");
}
