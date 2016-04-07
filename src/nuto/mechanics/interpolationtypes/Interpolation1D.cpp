/*
 * Interpolation1D.cpp
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/Interpolation1D.h"

NuTo::Interpolation1D::Interpolation1D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        InterpolationBase::InterpolationBase(rDofType, rTypeOrder, rDimension)
{

}

std::vector<Eigen::VectorXd> NuTo::Interpolation1D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    Eigen::VectorXd dummy; // has no influence in 1D
    return std::vector<Eigen::VectorXd>(1, this->CalculateNaturalSurfaceCoordinates(dummy, rSurface));
}

int NuTo::Interpolation1D::GetNumDofsPerNode() const
{
    switch (mDofType)
    {
    case NuTo::Node::COORDINATES:
        return mDimension;
    case NuTo::Node::DISPLACEMENTS:
        return mDimension;
    case NuTo::Node::TEMPERATURES:
        return 1;
    case NuTo::Node::NONLOCALEQSTRAIN:
        return 1;
    case NuTo::Node::NONLOCALEQPLASTICSTRAIN:
        return 2;
    case NuTo::Node::RELATIVEHUMIDITY:
        return 1;
    case NuTo::Node::WATERVOLUMEFRACTION:
        return 1;
    default:
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "dof type not found.");
    }
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
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InterpolationBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Interpolation1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation1D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Interpolation1D)
#endif  // ENABLE_SERIALIZATION
