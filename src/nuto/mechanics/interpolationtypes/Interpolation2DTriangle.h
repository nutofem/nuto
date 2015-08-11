/*
 * Interpolation2DTriangle.h
 *
 *  Created on: 20 Mar 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION2DTRIANGLE_H_
#define INTERPOLATION2DTRIANGLE_H_

#include "nuto/mechanics/interpolationtypes/Interpolation2D.h"

namespace NuTo
{
/**
@brief 2D triangular element with the following natural coordinate system and its surface parametrization
\f[\fbox{ \begin{tikzpicture}[scale=2]
  \draw[dotted, -latex] (0,0) -- (1.5,0) node[above]{$\xi$};
  \draw[dotted,-latex] (0,0) -- (0,1.5) node[above]{$\eta$};
  \draw[dashed] (0,0) node[below left]{$0$} -- (1,0) node[below]{$1$} -- (0,1) node[left]{$1$} -- (0,0);
  \draw[|-latex] (0.5,0) -- (0.8,0.0) node[midway, below] {$\alpha_0$};
  \draw[|-latex] (.5,.5) -- (0.2,0.8) node[midway, right] {$\alpha_1$};
  \draw[|-latex] (0,.5) -- (0.0,0.2) node[midway, left]  {$\alpha_2$};
  \node[rectangle, draw] at (.8,1) {$\alpha_i = [-1,1]$};
\end{tikzpicture}   } \f]
**/
class Interpolation2DTriangle: public Interpolation2D
{
public:
    Interpolation2DTriangle(const StructureBase* rStructure, NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    IntegrationType::eIntegrationType GetStandardIntegrationType() const override;

    //! @brief returns the natural coordinates of the dof node
    //! @param rDofType ... dof type
    //! @param rNodeIndexDof ... node index of the dof type
    const Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndexDof) const override;

    //! @brief calculates the shape functions for a specific dof
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    const Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    //! @brief returns derivative shape functions in the local coordinate system
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... map of derivative shape functions in the natural coordinate system for all dofs
    const Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    const Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    //! @brief returns the derivative of the surface parametrization
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... derivative of the surface parametrization
    const Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    //! @brief returns the number of surfaces
    int GetNumSurfaces() const override
    {
        return 3;
    }

private:
    //! @brief return the number node depending the shape and the order
    int CalculateNumNodes() const override;
};

} /* namespace NuTo */

#endif /* INTERPOLATIONTYPE2DTRIANGLE_H_ */
