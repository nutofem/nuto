// $Id$

#include "nuto/mechanics/cracks/CrackBase.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief Constructor
NuTo::CrackBase::CrackBase()
{}
NuTo::CrackBase::CrackBase(const std::string rName): mName(rName)
{}


//! @brief destructor
NuTo::CrackBase::~CrackBase()
{
}


#ifdef ENABLE_VISUALIZE
void NuTo::CrackBase::Visualize(VisualizeUnstructuredGrid& rVisualize) const
{
    throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] virtual function!");
}
#endif // ENABLE_VISUALIZE
