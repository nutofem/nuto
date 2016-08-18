//#ifdef ENABLE_SERIALIZATION
//#include <boost/archive/binary_oarchive.hpp>
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/xml_oarchive.hpp>
//#include <boost/archive/xml_iarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include "nuto/math/EigenBoostSerialization.h"
//#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/interpolationtypes/InterpolationBaseIGA.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

NuTo::InterpolationBaseIGA::InterpolationBaseIGA(Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder, int rDimension) :
    InterpolationBase::InterpolationBase(rDofType, rTypeOrder, rDimension)
{}


void NuTo::InterpolationBaseIGA::Initialize()
{
    mNumNodes = CalculateNumNodes();
    mNumDofs = mNumNodes*GetNumDofsPerNode();
}

#ifdef ENABLE_SERIALIZATION
NuTo::InterpolationBaseIGA::InterpolationBaseIGA():
    mTypeOrder(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1),
    mDofType(NuTo::Node::COORDINATES),
    mDimension(0)
{
}

template void NuTo::InterpolationBaseIGA::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::InterpolationBaseIGA::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::InterpolationBaseIGA::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::InterpolationBaseIGA::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::InterpolationBaseIGA::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::InterpolationBaseIGA::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::InterpolationBaseIGA::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize InterpolationBaseIGA" << std::endl;
#endif
    ar & boost::serialization::make_nvp("mTypeOrder", const_cast<NuTo::Interpolation::eTypeOrder&>(mTypeOrder));
    ar & boost::serialization::make_nvp("mDofType", const_cast<NuTo::Node::eDof&>(mDofType));
    ar & BOOST_SERIALIZATION_NVP(mIsConstitutiveInput);
    ar & BOOST_SERIALIZATION_NVP(mIsActive);
    ar & BOOST_SERIALIZATION_NVP(mNumDofs);
    ar & BOOST_SERIALIZATION_NVP(mNumNodes);

    ar & BOOST_SERIALIZATION_NVP(mLocalStartIndex);
    ar & boost::serialization::make_nvp("mDimension", const_cast<int&>(mDimension));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize InterpolationBaseIGA" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::InterpolationBaseIGA)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::InterpolationBaseIGA)
#endif  // ENABLE_SERIALIZATION
