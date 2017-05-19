/*
 * Interpolation2DQuad.h
 *
 *  Created on: 24 Apr 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/Interpolation2D.h"

namespace NuTo
{

/**
@brief 2D quadrilateral element with the following natural coordinate system and its surface parametrization
@image html Quad2D.png
**/
class Interpolation2DQuad: public Interpolation2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief just for serialization
    Interpolation2DQuad(){}
#endif  // ENABLE_SERIALIZATION
public:
    Interpolation2DQuad(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    eIntegrationType GetStandardIntegrationType() const override;

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    int GetNumSurfaces() const override
    {
        return 4;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:

    Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndexDof) const override;

    Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

private:
    
    int CalculateNumNodes() const override;

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Interpolation2DQuad)
#endif

