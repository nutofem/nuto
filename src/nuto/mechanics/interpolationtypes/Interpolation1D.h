/*
 * Interpolation1D.h
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"

namespace NuTo
{

class Interpolation1D: public InterpolationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief default constructor for serialization
protected:
    Interpolation1D(){}
#endif  // ENABLE_SERIALIZATION

public:

    Interpolation1D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    //! @brief return the number of dofs per node depending on dimension
    virtual int GetNumDofsPerNode() const override;

    //! @brief return the local dimension of the interpolation
    virtual int GetLocalDimension() const override
    {
        return 1;
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
BOOST_CLASS_EXPORT_KEY(NuTo::Interpolation1D)
#endif


