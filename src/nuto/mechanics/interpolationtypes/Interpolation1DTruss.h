/*
 * Interpolation1DTruss.h
 *
 *  Created on: 8 May 2015
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/interpolationtypes/Interpolation1D.h"

namespace NuTo
{

/**
@brief 2D quadrilateral element with the following natural coordinate system and its surface parametrization
\f[\fbox{ \begin{tikzpicture}
  \draw[dotted, |-latex] (0,0) -- (1.5,0) node[above]{$\xi$};
  \draw[dashed, |-|] (-1,0) node[below]{$-1$} node[above] {$\alpha_0$}
              --( 1,0) node[below]{$ 1$} node[above] {$\alpha_1$};
\end{tikzpicture}   } \f]
**/
class Interpolation1DTruss: public Interpolation1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief default constructor for serialization
protected:
    Interpolation1DTruss(){}
#endif  // ENABLE_SERIALIZATION

public:

    Interpolation1DTruss(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    IntegrationType::eIntegrationType GetStandardIntegrationType() const override;

    //! @brief returns the natural coordinates of the dof node
    //! @param rDofType ... dof type
    //! @param rNodeIndexDof ... node index of the dof type
    Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndexDof) const override;

    //! @brief calculates the shape functions for a specific dof
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    //! @brief returns derivative shape functions in the local coordinate system
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... map of derivative shape functions in the natural coordinate system for all dofs
    Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    //! @brief returns the derivative of the surface parametrization
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... derivative of the surface parametrization
    Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;


    //! @brief returns the number of surfaces
    inline int GetNumSurfaces() const override
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

private:

    //! @brief return the number node depending the shape and the order
    int CalculateNumNodes() const override;

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Interpolation1DTruss)
namespace boost
{
template<>
struct is_virtual_base_of<NuTo::Interpolation1D, NuTo::Interpolation1DTruss>: public mpl::true_ {};
}
#endif


