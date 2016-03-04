/*
 * Interpolation1D.h
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION1D_H_
#define INTERPOLATION1D_H_

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"

namespace NuTo
{

class Interpolation1D: public InterpolationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief default constructor for serialization
    Interpolation1D(){}

    Interpolation1D(const StructureBase* rStructure, NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    const std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    //! @brief return the number of dofs per node depending on dimension
    virtual int GetNumDofsPerNode() const override;

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

#endif /* INTERPOLATION1D_H_ */
