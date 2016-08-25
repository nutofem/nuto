/*
 * Interpolation3DBrick.h
 *
 *  Created on: 20 May 2015
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/interpolationtypes/Interpolation3D.h"

namespace NuTo
{
/**
@brief 3D tetrahedral element with the following natural coordinate system and its surface parametrization
\f[\fbox{
  \tikzset{engineering/.style={x={(-0.424cm,-0.265cm)}, y={(0.993cm,-0.122cm)}, z={(0cm,1cm)}}}
  \begin{tikzpicture}[engineering,scale = 1]
  \draw[-latex,dotted] ( 0.0 , 0.0 , 0.0 ) -- ( 1.0 , 0.0 , 0.0 ) node[pos=1.1] {$\xi$};
  \draw[-latex,dotted] ( 0.0 , 0.0 , 0.0 ) -- ( 0.0 , 1.0 , 0.0 ) node[pos=1.1] {$\eta$};
  \draw[-latex,dotted] ( 0.0 , 0.0 , 0.0 ) -- ( 0.0 , 0.0 , 1.0 ) node[pos=1.1] {$\zeta$};
  \draw[dashed] ( 1,-1,-1) -- (-1,-1,-1) -- (-1, 1,-1);
  \draw (-1, 1,-1) -- ( 1, 1,-1) -- ( 1,-1,-1);
  \draw (-1,-1, 1) -- ( 1,-1, 1) -- ( 1, 1, 1) -- (-1, 1, 1) -- cycle;
  \draw[dashed] (-1,-1, 1) -- (-1,-1,-1);
  \draw[] ( 1,-1, 1) -- ( 1,-1,-1);
  \draw[] ( 1, 1, 1) -- ( 1, 1,-1);
  \draw[] (-1, 1, 1) -- (-1, 1,-1);

  \begin{scope}[shift={(0,0,-2)}]
  \draw[dashed] (-1,-1,-1) -- ( 1,-1,-1) -- ( 1, 1,-1) -- (-1, 1,-1) -- cycle;
  \draw[-latex] (0,0,-1) -- (1,0,-1) node [pos=1.3] {$\beta_0$};
  \draw[-latex] (0,0,-1) -- (0,1,-1) node [pos=1.2] {$\alpha_0$};
  \end{scope}

  \begin{scope}[shift={(0,-2,0)}]
  \draw[dashed] (-1,-1,-1) -- ( 1,-1,-1) -- ( 1,-1, 1) -- (-1,-1, 1) -- cycle;
  \draw[-latex] (0,-1,0) -- (1,-1, 0) node [pos=1.3] {$\alpha_1$};
  \draw[-latex] (0,-1,0) -- (0,-1, 1) node [pos=1.2] {$\beta_1$};
  \end{scope}

  \begin{scope}[shift={(-2,0,0)}]
  \draw[dashed] (-1,-1,-1) -- (-1, 1,-1) -- (-1, 1, 1) -- (-1,-1, 1) -- cycle;
  \draw[-latex] (-1, 0,0) -- (-1, 1, 0) node [pos=1.2] {$\beta_2$};
  \draw[-latex] (-1,-0,0) -- (-1, 0, 1) node [pos=1.2] {$\alpha_2$};
  \end{scope}

    \begin{scope}[shift={(0,0,4)}]
  \draw[dashed] (-1,-1,-1) -- ( 1,-1,-1) -- ( 1, 1,-1) -- (-1, 1,-1) -- cycle;
  \draw[-latex] (0,0,-1) -- (1,0,-1) node [pos=1.3] {$\alpha_3$};
  \draw[-latex] (0,0,-1) -- (0,1,-1) node [pos=1.2] {$\beta_3$};
  \end{scope}

  \begin{scope}[shift={(0,4,0)}]
  \draw[dashed] (-1,-1,-1) -- ( 1,-1,-1) -- ( 1,-1, 1) -- (-1,-1, 1) -- cycle;
  \draw[-latex] (0,-1,0) -- (1,-1, 0) node [pos=1.3] {$\beta_4$};
  \draw[-latex] (0,-1,0) -- (0,-1, 1) node [pos=1.2] {$\alpha_4$};
  \end{scope}

  \begin{scope}[shift={(6,0,0)}]
  \draw[dashed] (-1,-1,-1) -- (-1, 1,-1) -- (-1, 1, 1) -- (-1,-1, 1) -- cycle;
  \draw[-latex] (-1, 0,0) -- (-1, 1, 0) node [pos=1.2] {$\alpha_5$};
  \draw[-latex] (-1,-0,0) -- (-1, 0, 1) node [pos=1.2] {$\beta_5$};
  \end{scope}

  \node at (-2,1,-2.5) {$\mathbf{ \alpha} \times \mathbf{\beta}$ pointing outwards};
\end{tikzpicture}   } \f]
**/
class Interpolation3DBrick: public Interpolation3D
{
public:
    Interpolation3DBrick(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

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

    //! @brief returns the number of surfaces
    int GetNumSurfaces() const override
    {
        return 6;
    }

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

private:
    //! @brief return the number node depending the shape and the order
    int CalculateNumNodes() const override;
};

} /* namespace NuTo */

