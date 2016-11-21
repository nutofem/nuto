/*
 * InterpolationBaseFEM.h
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/MechanicsException.h"

#include <vector>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "nuto/math/CustomBoostSerializationExtensions.h"
#endif  // ENABLE_SERIALIZATION

namespace NuTo
{
class IntegrationTypeBase;

//! @brief this class stores the information of the interpolation of a single dof type
//! @remark the API only allows const access to this class via the InterpolationType.Get(dofType)
//! method. Its data members are set via the friend class property.
class InterpolationBaseFEM : public InterpolationBase
{
friend class InterpolationType;

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    InterpolationBaseFEM();
#endif

public:
    InterpolationBaseFEM(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    virtual ~InterpolationBaseFEM() {}

    virtual int GetSplineDegree(int dir) const override;

    //********************************************
    //             NODE METHODS
    //********************************************

    //! @brief returns the natural coordinates of the dof node
    //! @param rNodeIndex ... node index
    const Eigen::VectorXd& GetNaturalNodeCoordinates(int rNodeIndex) const;

    //! @brief returns the natural coordinates of the dof node
    //! @param rNodeIndex ... node index
    virtual Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndex) const = 0;

    void CalculateSurfaceNodeIds();

    //********************************************
    //       SHAPE FUNCTIONS
    //********************************************

    //! @brief returns specific shape functions via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific shape functions
    const Eigen::VectorXd& GetShapeFunctions(int rIP) const;

    //! @brief returns specific N-matrix via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific N-matrix
    const Eigen::MatrixXd& GetMatrixN(int rIP) const;

    //! @brief calculates the shape functions for a specific dof
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    virtual Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const = 0;

    //! @brief calculates the N-Matrix, blows up the shape functions to the correct format (e.g. 3D: N & 0 & 0 \\ 0 & N & 0 \\ 0 & 0 & N ...)
    Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const;


    // --- IGA interpolation--- //

    //! @brief returns specific shape functions at a parameter, whicg fits the knot vector
    //! @param rIP ... id of the integration point
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (a transformation needs to be done, since integration point coordinates are in [-1, 1])
    //! @return ... specific shape functions
    virtual Eigen::VectorXd CalculateShapeFunctions(int rIP, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    //! @brief returns the N matrix at a parameter for IGA elements
    //! @param rIP ... id of the integration point
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateMatrixN(int rIP, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    //! @brief returns the N matrix for IGA elements at a parameter, which fits to the knot vector (e.g. 3D: N & 0 & 0 \\ 0 & N & 0 \\ 0 & 0 & N ...)
    //! @param rCoordinates ... parameter
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    //! @brief returns the N matrix at a parameter, which fits to the knot vector (e.g. 3D: N & 0 & 0 \\ 0 & N & 0 \\ 0 & 0 & N ...)
    //! @param rParameters ... parameter on the curve
    //! @param rKnotIDs ... knot span
    //! @param rDerivative ... the order of derivative (only 0,1,2 possible)
    //! @param rDirection ... for 1D only 0 (in 2D 0(x) and 1(y))
    virtual Eigen::MatrixXd CalculateMatrixNDerivative(const Eigen::VectorXd& rParameters, const Eigen::VectorXi& rKnotIDs, int rDerivative, int rDirection) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "So far implemeneted only for IGA!");
    }
    //********************************************
    //       DERIVATIVE SHAPE FUNCTIONS NATURAL
    //********************************************

    //! @brief returns specific derivative shape functions natural via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific derivative shape functions natural
    const Eigen::MatrixXd & GetDerivativeShapeFunctionsNatural(int rIP) const override;

    //! @brief returns specific derivative shape functions natural via coordinates
    //! @param rCoordinates ... integration point coordinates
    //! @return ... specific derivative shape functions natural
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override = 0;

    // --- IGA interpolation--- //

    //! @brief returns specific derivative shape functions at a parameter, which fits to the knot vector
    //! @param rIP ... id of the integration point
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(int rIP, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    //! @brief returns specific derivative shape functions at a parameter, which fits to the knot vector
    //! @param rCoordinates ... parameter
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    //********************************************
    //       SURFACE PARAMETRIZATION
    //********************************************

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const = 0;

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface, const Eigen::MatrixXd &rKnots) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

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

    Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    int GetSurfaceDegree(int rSurface) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class, this is the load routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    //! @brief serializes the class, this is the save routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

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
    bool NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const;

    //! @brief calculate and store the shape functions and their derivatives
    //! @param rIntegrationType ... integration type
    void UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType) override;

    //! @brief return the number node depending the shape and the order
    virtual int CalculateNumNodes() const = 0;

    //! @brief this method sets the mNumDofs, mNumNodes and mNodeIndices members
    //! @remark it should be called from the ctor InterpolationTypeBase()
    //! but uses pure virutal functions. Thus it must be called in the
    //! ctors of the child classes.
    void Initialize();

    //********************************************
    //               MEMBERS
    //********************************************

    // members for each integration point
    std::vector<Eigen::VectorXd> mNodeCoordinates;
    std::vector<Eigen::VectorXd> mShapeFunctions;
    std::vector<Eigen::MatrixXd> mMatrixN;
    std::vector<Eigen::MatrixXd> mDerivativeShapeFunctionsNatural;
};
} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::InterpolationBaseFEM)
#endif
