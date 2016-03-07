/*
 * InterpolationType2D.cpp
 *
 *  Created on: 23 Mar 2015
 *      Author: ttitsche
 */

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/interpolationtypes/Interpolation2D.h"

NuTo::Interpolation2D::Interpolation2D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        InterpolationBase::InterpolationBase(rDofType, rTypeOrder, rDimension)
{

}

const std::vector<Eigen::VectorXd> NuTo::Interpolation2D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    int numNodes = 2;
    // returns exactly two nodes, one at alpha = -1 and one at alpha = 1 of the required surface
    Eigen::Vector2d alpha(-1, 1);
    std::vector<Eigen::VectorXd> surfaceEdgeCoordinates(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
        Eigen::VectorXd naturalSurfaceCoordinate = alpha.row(i);
        surfaceEdgeCoordinates[i] = CalculateNaturalSurfaceCoordinates(naturalSurfaceCoordinate, rSurface);
    }
    return surfaceEdgeCoordinates;
}

int NuTo::Interpolation2D::GetNumDofsPerNode() const
{
    switch (mDofType)
    {
    case NuTo::Node::COORDINATES:
        return 2;
    case NuTo::Node::DISPLACEMENTS:
        return 2;
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
        throw NuTo::MechanicsException("[NuTo::Interpolation2D::GetNumDofsPerNode] dof type not found.");
    }
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::Interpolation2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Interpolation2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Interpolation2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Interpolation2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Interpolation2D\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InterpolationBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Interpolation2D\n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation2D)
#endif
