// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constraints/ConstraintBase.h"

//! @brief constructor
NuTo::ConstraintBase::ConstraintBase()
{
}

// destructor
NuTo::ConstraintBase::~ConstraintBase()
{
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintBase::serialize(Archive & ar, const unsigned int version)
{
}
#endif // ENABLE_SERIALIZATION

//!@brief sets/modifies the right hand side of the constraint equations
//!@param rRHS new right hand side
void NuTo::ConstraintBase::SetRHS(double rRHS)
{
    throw MechanicsException("[NuTo::ConstraintBase] Set right hand side for this type of constraints not implemented.");
}

//!@brief set the strain of the periodic boundary conditions
//!@param rStrain strain (e_xx,e_yy,gamma_xy)
void NuTo::ConstraintBase::SetStrain(const NuTo::FullMatrix<double>& rStrain)
{
    throw MechanicsException("[NuTo::ConstraintBase] Set strain for this type of constraints not implemented.");
}
