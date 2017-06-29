/*
 * Interpolation1D.cpp
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#include "base/Exception.h"
#include "mechanics/interpolationtypes/Interpolation1D.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

NuTo::Interpolation1D::Interpolation1D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        InterpolationBaseFEM::InterpolationBaseFEM(rDofType, rTypeOrder, rDimension)
{

}

std::vector<Eigen::VectorXd> NuTo::Interpolation1D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    Eigen::VectorXd dummy; // has no influence in 1D
    return std::vector<Eigen::VectorXd>(1, this->CalculateNaturalSurfaceCoordinates(dummy, rSurface));
}


#ifdef ENABLE_SERIALIZATION
template void NuTo::Interpolation1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Interpolation1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Interpolation1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InterpolationBaseFEM);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Interpolation1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation1D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Interpolation1D)
#endif  // ENABLE_SERIALIZATION
