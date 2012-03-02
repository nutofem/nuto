// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#include <boost/ptr_container/serialize_ptr_list.hpp>
#endif // ENABLE_SERIALIZATION

#ifdef SHOW_TIME
#include <ctime>
#endif

# ifdef _OPENMP
#include <omp.h>
# endif


#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/mechanics/MechanicsException.h"


NuTo::TimeIntegrationBase::TimeIntegrationBase()  : NuTo::NuToObject::NuToObject()
{
	mConstraintLoad = -1;
	mConstraintRHSDelta = 0.;
	mTime = 0.;
	mTimeDelta = 1.;
	mLoadStep = 1;
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::TimeIntegrationBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::TimeIntegrationBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::TimeIntegrationBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    mLogger << "start serialization of TimeIntegrationBase" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP(mConstraintLoad)
       & BOOST_SERIALIZATION_NVP(mConstraintRHSDelta)
       & BOOST_SERIALIZATION_NVP(mTime)
       & BOOST_SERIALIZATION_NVP(mTimeDelta)
       & BOOST_SERIALIZATION_NVP(mLoadStep);
#ifdef DEBUG_SERIALIZATION
    mLogger << "finish serialization of structure base" << "\n";
#endif
}
#endif  // ENABLE_SERIALIZATION

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::TimeIntegrationBase::Info()const
{
}

