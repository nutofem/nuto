/*
 * DofStatus.cpp
 *
 *  Created on: 4 Apr 2016
 *      Author: ttitsche
 */

#include "nuto/mechanics/dofSubMatrixStorage/DofStatus.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/map.hpp>
#endif // ENABLE_SERIALIZATION


NuTo::DofStatus::DofStatus() : mHasInteractingConstraints(false)
{
}

namespace NuTo
{
std::ostream& operator<<(std::ostream& out, const DofStatus& dofStatus)
{
    out << "[---DofStatus:---]\n";
    out << "Existing DOF types:\n";
    for(auto dof : dofStatus.mDofTypes)
    {
        out << dof << "; ";
    }
    out << "Active DOF types:\n";
    for(auto dof : dofStatus.mActiveDofTypes)
    {
        out << dof << "; ";
    }
    for(auto dof : dofStatus.mNumActiveDofs)
    {
        out << "Number of active Dofs of type " << dof.first << ": " << dof.second << "\n";
    }
    for(auto dof : dofStatus.mNumDependentDofs)
    {
        out << "Number of dependent Dofs of type " << dof.first << ": " << dof.second << "\n";
    }
    return out;
}
} // namespace NuTo

#ifdef ENABLE_SERIALIZATION
template void NuTo::DofStatus::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::DofStatus::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::DofStatus::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::DofStatus::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::DofStatus::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::DofStatus::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::DofStatus::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of DofStatus" << "\n";
#endif
    ar & BOOST_SERIALIZATION_NVP(mNumActiveDofs)
       & BOOST_SERIALIZATION_NVP(mNumDependentDofs)
       & BOOST_SERIALIZATION_NVP(mDofTypes)
       & BOOST_SERIALIZATION_NVP(mActiveDofTypes)
       & BOOST_SERIALIZATION_NVP(mSymmetricDofTypes)
       & BOOST_SERIALIZATION_NVP(mHasInteractingConstraints);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of DofStatus" << "\n";
#endif
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::DofStatus)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
