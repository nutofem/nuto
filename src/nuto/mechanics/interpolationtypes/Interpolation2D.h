/*
 * Interpolation2D.h
 *
 *  Created on: 23 Mar 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION2D_H_
#define INTERPOLATION2D_H_

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"

namespace NuTo
{

class Interpolation2D: public InterpolationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief just for serialization
    Interpolation2D(){}

    Interpolation2D(const StructureBase* rStructure, NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    const std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    //! @brief return the number of dofs per node depending on dimension
    int GetNumDofsPerNode() const override;

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

#endif /*INTERPOLATION2D_H_ */
