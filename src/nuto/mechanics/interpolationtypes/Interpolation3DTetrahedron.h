/*
 * Interpolation3DTetrahedron.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#ifndef INTERPOLATION3DTETRAHEDRON_H_
#define INTERPOLATION3DTETRAHEDRON_H_

#include "nuto/mechanics/interpolationtypes/Interpolation3D.h"

namespace NuTo
{
/**
@brief 3D tetrahedral element with the following natural coordinate system and its surface parametrization
\f[\fbox{
  \tikzset{engineering/.style={x={(-0.424cm,-0.265cm)}, y={(0.993cm,-0.122cm)}, z={(0cm,1cm)}}}
  \begin{tikzpicture}[engineering,scale = 2]
  \draw[-latex,dotted] ( 0.0 , 0.0 , 0.0 ) -- ( 1.5 , 0.0 , 0.0 ) node[pos=1.1] {$\xi$};
  \draw[-latex,dotted] ( 0.0 , 0.0 , 0.0 ) -- ( 0.0 , 1.5 , 0.0 ) node[pos=1.1] {$\eta$};
  \draw[-latex,dotted] ( 0.0 , 0.0 , 0.0 ) -- ( 0.0 , 0.0 , 1.5 ) node[pos=1.1] {$\zeta$};
  \draw[dashed] (0,0,0) -- (1,0,0);
  \draw[dashed] (0,0,0) -- (0,1,0);
  \draw[dashed] (0,0,0) -- (0,0,1);
  \draw[] (1,0,0) -- (0,1,0) -- (0,0,1) -- cycle;

  \begin{scope}[shift={(0,0,-.7)}]
      \draw[dashed] (0,0,0) -- (1,0,0) -- (0,1,0) -- cycle;
      \draw[-latex] (0,0,0) -- (1,0,0) node [pos=1.2] {$\beta_0$};
      \draw[-latex] (0,0,0) -- (0,1,0) node [pos=1.15] {$\alpha_0$};
  \end{scope}

  \begin{scope}[shift={(-3,0,0)}]
      \draw[dashed] (0,0,0) -- (0,1,0) -- (0,0,1) -- cycle;
      \draw[-latex] (0,0,0) -- (0,1,0) node [pos=1.15] {$\beta_1$};
      \draw[-latex] (0,0,0) -- (0,0,1) node [pos=1.15] {$\alpha_1$};
  \end{scope}

  \begin{scope}[shift={(0,-.7,0)}]
      \draw[dashed] (0,0,0) -- (1,0,0) -- (0,0,1) -- cycle;
      \draw[-latex] (0,0,0) -- (1,0,0) node [pos=1.2] {$\alpha_2$};
      \draw[-latex] (0,0,0) -- (0,0,1) node [pos=1.15] {$\beta_2$};
  \end{scope}

  \begin{scope}[shift={(1,1,1)}]
      \draw[dashed] (1,0,0) -- (0,1,0) -- (0,0,1) -- cycle;
      \draw[-latex] (1,0,0) -- (0,1,0) node [pos=1.1] {$\alpha_3$};
      \draw[-latex] (1,0,0) -- (0,0,1) node [pos=1.1] {$\beta_3$};
  \end{scope}

  \node at (-2,1,-1) {$\mathbf{ \alpha} \times \mathbf{\beta}$ pointing outwards};
\end{tikzpicture}   } \f]
**/
class Interpolation3DTetrahedron: public Interpolation3D
{
public:
    Interpolation3DTetrahedron(const StructureBase* rStructure, NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder);

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
        return 4;
    }

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    const std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

private:
    //! @brief return the number node depending the shape and the order
    int CalculateNumNodes() const override;
};

} /* namespace NuTo */

#endif /* INTERPOLATION3DTETRAHEDRON_H_ */
