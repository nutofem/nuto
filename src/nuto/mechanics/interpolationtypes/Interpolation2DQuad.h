/*
 * Interpolation2DQuad.h
 *
 *  Created on: 24 Apr 2015
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/interpolationtypes/Interpolation2D.h"

namespace NuTo
{

/**
@brief 2D quadrilateral element with the following natural coordinate system and its surface parametrization
\f[\fbox{ \begin{tikzpicture}
  \draw[dotted, -latex] (0,0) -- (1.5,0) node[above]{$\xi$};
  \draw[dotted,-latex] (0,0) -- (0,1.5) node[above]{$\eta$};
  \draw[dashed] (-1,-1) node[below  left]{$(-1,-1)$}
              --( 1,-1) node[below right]{$( 1,-1)$}
              --( 1, 1) node[above right]{$( 1, 1)$}
              --(-1, 1) node[above  left]{$(-1, 1)$} -- cycle;
  \draw[|-latex] ( 0,-1) -- (.8,-1) node[midway, below] {$\alpha_0$};
  \draw[|-latex] ( 1, 0) -- ( 1,.8) node[midway, above right] {$\alpha_1$};
  \draw[|-latex] ( 0, 1) -- (-.8, 1) node[midway, below] {$\alpha_2$};
  \draw[|-latex] (-1, 0) -- (-1,-.8) node[midway, left] {$\alpha_3$};
  \node[rectangle, draw] at (1.5,2) {$\alpha_i = [-1,1]$};
\end{tikzpicture}   } \f]
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

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    eIntegrationType GetStandardIntegrationType() const override;

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

    //! @brief updates the mNodeIndices member according to the order of interpolation
    void UpdateNodeIndices(const std::vector<Eigen::VectorXd> &rNodeCoordinates, std::function<bool(const Eigen::VectorXd& rC1, const Eigen::VectorXd& rC2)> rFunction) override;

    //! @brief returns the number of surfaces
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

private:
    //! @brief return the number node depending the shape and the order
    int CalculateNumNodes() const override;

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Interpolation2DQuad)
#endif

