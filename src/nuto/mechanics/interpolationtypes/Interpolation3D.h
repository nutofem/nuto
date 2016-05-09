/*
 * Interpolation3D.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"

namespace NuTo
{

class Interpolation3D: public InterpolationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief just for serialization
protected:
    Interpolation3D(){}
#endif  // ENABLE_SERIALIZATION
public:
    Interpolation3D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Interpolation3D" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InterpolationBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Interpolation3D" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

    //! @brief return the number of dofs per node depending on dimension
    virtual int GetNumDofsPerNode() const override;

    //! @brief return the local dimension of the interpolation
    virtual int GetLocalDimension() const override
    {
        return 3;
    }

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Interpolation3D)
#endif
