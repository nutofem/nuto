/*
 * Interpolation2D.h
 *
 *  Created on: 23 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/InterpolationBaseFEM.h"

namespace NuTo
{

class Interpolation2D: public InterpolationBaseFEM
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief just for serialization
protected:
    Interpolation2D(){}
#endif  // ENABLE_SERIALIZATION

public:
    Interpolation2D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const override
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Implemented in NuTo::InterpolationType::GetSurfaceNodeIndices.");
    }

    virtual int GetLocalDimension() const override
    {
        return 2;
    }
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Interpolation2D)
#endif

