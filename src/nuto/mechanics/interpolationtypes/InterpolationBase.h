/*
 * InterpolationBase.h
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/mechanics/structures/StructureBase.h"


#include <vector>
#include <assert.h>
#include <string>
#include <sstream>

namespace NuTo
{
class IntegrationTypeBase;

//! @brief this class stores the information of the interpolation of a single dof type
//! @remark the API only allows const access to this class via the InterpolationType.Get(dofType)
//! method. Its data members are set via the friend class property.
class InterpolationBase
{
friend class InterpolationType;

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    InterpolationBase();
#endif

public:
    InterpolationBase(NuTo::Node::eDof rDofType, Interpolation::eTypeOrder rTypeOrder, int rDimension);

    virtual ~InterpolationBase() {}

    virtual NuTo::Interpolation::eTypeOrder GetTypeOrder() const;

    virtual int GetSplineDegree() const = 0;

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    virtual IntegrationType::eIntegrationType GetStandardIntegrationType() const = 0;

    //********************************************
    //             DOF METHODS
    //********************************************

    //! @brief returns whether or not the dof is active
    bool IsActive() const;

    //! @brief returns whether or not the dof is constitutive input
    bool IsConstitutiveInput() const;

    //! @brief returns the number of dofs
    int GetNumDofs() const;

    //! @brief returns the start index that describes where to place the element sub-vectors and sub-matrices of rDofType
    //! in a multi-physics context in the element vectors and matrices
    int GetLocalStartIndex() const;

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
    virtual const Eigen::VectorXd& GetShapeFunctions(int rIP) const = 0;

    //! @brief returns specific N-matrix via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific N-matrix
    virtual const Eigen::MatrixXd& GetMatrixN(int rIP) const = 0;

    //! @brief calculates the shape functions for a specific dof
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    virtual Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const = 0;

    //! @brief calculates the N-Matrix, blows up the shape functions to the correct format (e.g. 3D: N & 0 & 0 \\ 0 & N & 0 \\ 0 & 0 & N ...)
    virtual Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const = 0;

    //! @brief calculates the shape functions for a specific dof
    //! @param rIP ... integration point index
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    virtual Eigen::VectorXd CalculateShapeFunctions(int rIP, const Eigen::VectorXi &mKnotIndices) const = 0;

    //! @brief calculates the N-Matrix, blows up the shape functions to the correct format (e.g. 3D: N & 0 & 0 \\ 0 & N & 0 \\ 0 & 0 & N ...)
    virtual Eigen::MatrixXd CalculateMatrixN(int rIP, const Eigen::VectorXi &mKnotIndices) const = 0;

    //********************************************
    //       DERIVATIVE SHAPE FUNCTIONS NATURAL
    //********************************************

    //! @brief returns specific derivative shape functions natural via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific derivative shape functions natural
    virtual const Eigen::MatrixXd & GetDerivativeShapeFunctionsNatural(int rIP) const = 0;

    //! @brief returns specific derivative shape functions natural via coordinates
    //! @param rCoordinates ... integration point coordinates
    //! @return ... specific derivative shape functions natural
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const = 0;

    //! @brief returns specific derivative shape functions natural via coordinates
    //! @param rIP ... integration point id
    //! @return ... specific derivative shape functions natural
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(int rIP, const Eigen::VectorXi &mKnotIndices) const = 0;

    //********************************************
    //       SURFACE PARAMETRIZATION
    //********************************************

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const = 0;

    //! @brief returns the derivative of the surface parametrization
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... derivative of the surface parametrization
    virtual Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const = 0;

    //! @brief returns the number of surfaces
    virtual int GetNumSurfaces() const = 0;

    //! @brief return the number of dofs per node depending on dimension
    virtual int GetNumDofsPerNode() const = 0;

    //! @brief return the local dimension of the interpolation
    virtual int GetLocalDimension() const = 0;


#ifdef ENABLE_SERIALIZATION
//    //! @brief serializes the class, this is the load routine
//    //! @param ar         archive
//    //! @param version    version
//    template<class Archive>
//    void load(Archive & ar, const unsigned int version);

//    //! @brief serializes the class, this is the save routine
//    //! @param ar         archive
//    //! @param version    version
//    template<class Archive>
//    void save(Archive & ar, const unsigned int version) const;

//    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

#endif  // ENABLE_SERIALIZATION

protected:

    //! @brief returns the natural coordinates of the nodes that span the surface
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural surface edge coordinates
    virtual std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const = 0;

    //! @brief returns true if a node is on the surface
    //! @param rSurface ... surface id
    //! @param rNaturalNodeCoordinate ... natural coordinate of the node to test
    virtual bool NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const = 0;

    //! @brief calculate and store the shape functions and their derivatives
    //! @param rIntegrationType ... integration type
    virtual void UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType) {}

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

    int mLocalStartIndex;

    // members for each surface
    std::vector<std::vector<int>> mSurfaceNodeIndices;

    bool mUpdateRequired;

    //! @brief dimension = Structure.GetDimension()
    const int mDimension;
};
} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::InterpolationBase)
#endif


