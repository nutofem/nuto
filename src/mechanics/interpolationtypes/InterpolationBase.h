/*
 * InterpolationBase.h
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include <vector>
#include <Eigen/Core>

namespace NuTo
{
enum class eIntegrationType;
namespace Interpolation
{
enum class eTypeOrder;
} // namespace Interpolation

namespace Node
{
enum class eDof : unsigned char;
} // namespace Node

//! @brief this class stores the information of the interpolation of a single dof type
//! @remark the API only allows const access to this class via the InterpolationType.Get(dofType)
//! method. Its data members are set via the friend class property.
class InterpolationBase
{
    friend class InterpolationType;

public:
    InterpolationBase(NuTo::Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder, int rDimension);

    virtual ~InterpolationBase()
    {
    }

    virtual NuTo::Interpolation::eTypeOrder GetTypeOrder() const;

    virtual int GetSplineDegree(int dir) const = 0;

    virtual void ClearCache() const
    {
    }

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    virtual eIntegrationType GetStandardIntegrationType() const = 0;

    //********************************************
    //             DOF METHODS
    //********************************************

    //! @brief returns whether or not the dof is active
    bool IsActive() const;

    //! @brief returns whether or not the dof is constitutive input
    bool IsConstitutiveInput() const;

    //! @brief returns the number of dofs
    int GetNumDofs() const;

    //********************************************
    //             NODE METHODS
    //********************************************

    //! @brief returns the total number of nodes
    int GetNumNodes() const;

    //! @brief returns the node index of a specific DOF node
    //! @param rNodeDofIndex ... node dof index
    int GetNodeIndex(int rNodeDofIndex) const;

    //! @brief returns the total number of nodes on a surface
    //! @param rSurface ... surface id
    int GetNumSurfaceNodes(int rSurface) const;

    //! @brief returns the node index of a specific DOF node on the surface
    //! @param rSurface ... surface id
    //! @param rNodeDofIndex ... node dof index
    int GetSurfaceNodeIndex(int rSurface, int rNodeDofIndex) const;

    //! @brief returns the natural coordinates of the dof node
    //! @param rNodeIndex ... node index
    virtual const Eigen::VectorXd& GetNaturalNodeCoordinates(int rNodeIndex) const = 0;

    //! @brief returns the natural coordinates of the dof node
    //! @param rNodeIndex ... node index
    virtual Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndex) const = 0;

    virtual void CalculateSurfaceNodeIds() = 0;

    //********************************************
    //       SHAPE FUNCTIONS
    //********************************************

    //! @brief returns specific shape functions via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific shape functions
    virtual const Eigen::VectorXd& ShapeFunctions(const Eigen::VectorXd& naturalCoordinates) const = 0;

    //! @brief returns specific N-matrix via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific N-matrix
    virtual const Eigen::MatrixXd& MatrixN(const Eigen::VectorXd& naturalCoordinates) const = 0;


    // --- IGA interpolation--- //

    //! @brief returns specific shape functions at a parameter, whicg fits the knot vector
    //! @param rIP ... id of the integration point
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (a transformation needs
    //! to be done, since integration point coordinates are in [-1, 1])
    //! @return ... specific shape functions
    virtual Eigen::VectorXd ShapeFunctionsIGA(const Eigen::VectorXd& naturalCoordinates,
                                              const Eigen::VectorXi& rKnotIDs) const = 0;

    //! @brief returns the N matrix for IGA elements at a parameter, which fits to the knot vector (e.g. 3D: N & 0 & 0
    //! \\ 0 & N & 0 \\ 0 & 0 & N ...)
    //! @param rCoordinates ... parameter
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd MatrixNIGA(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi& rKnotIDs) const = 0;

    //! @brief returns the N matrix at a parameter, which fits to the knot vector (e.g. 3D: N & 0 & 0 \\ 0 & N & 0 \\ 0
    //! & 0 & N ...)
    //! @param rParameters ... parameter on the curve
    //! @param rKnotIDs ... knot span
    //! @param rDerivative ... the order of derivative (only 0,1,2 possible)
    //! @param rDirection ... for 1D only 0 (in 2D 0(x) and 1(y))
    virtual Eigen::MatrixXd MatrixNDerivativeIGA(const Eigen::VectorXd& rParameters, const Eigen::VectorXi& rKnotIDs,
                                                 int rDerivative, int rDirection) const = 0;

    //********************************************
    //       DERIVATIVE SHAPE FUNCTIONS NATURAL
    //********************************************

    //! @brief returns specific derivative shape functions natural via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific derivative shape functions natural
    virtual const Eigen::MatrixXd& DerivativeShapeFunctionsNatural(const Eigen::VectorXd& naturalCoordinates) const = 0;

    // --- IGA interpolation--- //

    //! @brief returns specific derivative shape functions at a parameter, which fits to the knot vector
    //! @param rCoordinates ... parameter
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd DerivativeShapeFunctionsNaturalIGA(const Eigen::VectorXd& rCoordinates,
                                                               const Eigen::VectorXi& rKnotIDs) const = 0;

    //********************************************
    //       SURFACE PARAMETRIZATION
    //********************************************

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                               int rSurface) const = 0;

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @param rKnots ... knots bounding the IGA element
    //! @return ... natural coordinates of the elements surface
    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinatesIGA(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                                  int rSurface,
                                                                  const Eigen::MatrixXd& rKnots) const = 0;

    //! @brief returns the derivative of the surface parametrization
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... derivative of the surface parametrization
    virtual Eigen::MatrixXd
    CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                 int rSurface) const = 0;

    //! @brief returns the number of surfaces
    virtual int GetNumSurfaces() const = 0;

    //! @brief return the local dimension of the interpolation
    virtual int GetLocalDimension() const = 0;

    virtual Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const = 0;

    virtual int GetSurfaceDegree(int rSurface) const = 0;

protected:
    //! @brief calculates the shape functions for a specific dof
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    virtual Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const = 0;

    //! @brief calculates the N-Matrix, blows up the shape functions to the correct format (e.g. 3D: N & 0 & 0 \\ 0 & N
    //! & 0 \\ 0 & 0 & N ...)
    virtual Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const = 0;

    //! @brief returns specific derivative shape functions natural via coordinates
    //! @param rCoordinates ... integration point coordinates
    //! @return ... specific derivative shape functions natural
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const = 0;

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    virtual std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const = 0;

    //! @brief returns true if a node is on the surface
    //! @param rSurface ... surface id
    //! @param rNaturalNodeCoordinate ... natural coordinate of the node to test
    virtual bool NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const = 0;

    //! @brief return the number node depending the shape and the order
    virtual int CalculateNumNodes() const = 0;

    //********************************************
    //               MEMBERS
    //********************************************

    // dof members - simple storage
    const NuTo::Node::eDof mDofType;

    bool mIsConstitutiveInput;
    bool mIsActive;

    int mNumDofs;
    int mNumNodes;

    const NuTo::Interpolation::eTypeOrder mTypeOrder;

    std::vector<int> mNodeIndices;

    // members for each surface
    std::vector<std::vector<int>> mSurfaceNodeIndices;

    //! @brief dimension = Structure.GetDimension()
    const int mDimension;
};
} /* namespace NuTo */
