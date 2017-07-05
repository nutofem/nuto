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

#include "base/Exception.h"
#include "mechanics/interpolationtypes/Interpolation2D.h"

NuTo::Interpolation2D::Interpolation2D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder,
                                       int rDimension)
    : InterpolationBaseFEM::InterpolationBaseFEM(rDofType, rTypeOrder, rDimension)
{
}

std::vector<Eigen::VectorXd> NuTo::Interpolation2D::GetSurfaceEdgesCoordinates(int rSurface) const
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
