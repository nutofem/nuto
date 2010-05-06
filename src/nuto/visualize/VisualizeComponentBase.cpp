// $ld: $ 
// VisualizeComponentBase.cpp
// created Apr 27, 2010 by Joerg F. Unger

#include "nuto/mechanics/elements/ElementWithDataBase.h"
#include "nuto/visualize/VisualizeComponentBase.h"
#include "nuto/visualize/VisualizeException.h"

int NuTo::VisualizeComponentBase::GetElementId()const
{
	throw VisualizeException("[NuTo::VisualizeComponentBase::GetElementId] Visualization component has no ElementId.");
}

const NuTo::ElementWithDataBase* NuTo::VisualizeComponentBase::GetElement()const
{
	throw VisualizeException("[NuTo::VisualizeComponentBase::GetElement] Visualization component has no Element.");
}

int NuTo::VisualizeComponentBase::GetIp()const
{
	throw VisualizeException("[NuTo::VisualizeComponentBase::GetIp] Visualization component has no Integration point.");
}


