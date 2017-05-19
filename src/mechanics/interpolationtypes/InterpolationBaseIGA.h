#pragma once

#include "mechanics/interpolationtypes/InterpolationBase.h"

#include "mechanics/MechanicsException.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "math/CustomBoostSerializationExtensions.h"
#endif  // ENABLE_SERIALIZATION

#include <vector>

namespace NuTo
{
class IntegrationTypeBase;

//! @brief this class stores the information of the interpolation of a single dof type
//! @remark the API only allows const access to this class via the InterpolationType.Get(dofType)
//! method. Its data members are set via the friend class property.
class InterpolationBaseIGA : public InterpolationBase
{
friend class InterpolationType;

#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
protected:
    InterpolationBaseIGA();
#endif

public:
    InterpolationBaseIGA(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    virtual ~InterpolationBaseIGA() {}

    virtual int GetSplineDegree(int dir) const override = 0;

    //! @brief determines the standard integration type depending on shape, type and order
    //! @return standard integration type
    virtual eIntegrationType GetStandardIntegrationType() const override = 0;

    //********************************************
    //             NODE METHODS
    //********************************************

    //! @brief returns the natural coordinates of the dof node
    //! @param rNodeIndex ... node index
    const Eigen::VectorXd& GetNaturalNodeCoordinates(int rNodeIndex) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "No natural node coordinates in IGA!");
    }

    //! @brief returns the natural coordinates of the dof node
    //! @param rNodeIndex ... node index
    virtual Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndex) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "No natural node coordinates in IGA!");
    }

    void CalculateSurfaceNodeIds() override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "No natural node coordinates in IGA!");
    }

    //********************************************
    //       SHAPE FUNCTIONS
    //********************************************

    //! @brief returns specific shape functions via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific shape functions
    virtual const Eigen::VectorXd& GetShapeFunctions(int rIP) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "The shape functions are calculated on the fly, use 'GetShapeFunctions' routine!");
    }

    //! @brief returns specific N-matrix via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific N-matrix
    virtual const Eigen::MatrixXd& GetMatrixN(int rIP) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "The shape functions are calculated on the fly, use 'GetMatrixN' routine!");
    }

    // -- shape functions --//

    //! @brief calculates the shape functions for a specific dof
    //! @param rCoordinates ... integration point coordinates
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    virtual Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override = 0;

    //! @brief calculates the shape functions for a specific dof
    //! @param rIP ... integration point index
    //! @param rDofType ... dof type
    //! @return ... shape functions for the specific dof type
    virtual Eigen::VectorXd CalculateShapeFunctions(int rIP, const Eigen::VectorXi &rKnotIDs) const override = 0;

    // -- N matrix --//

    //! @brief returns the N matrix at a parameter for IGA elements
    //! @param rIP ... id of the integration point
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateMatrixN(int rIP, const Eigen::VectorXi &rKnotIDs) const override= 0;

    //! @brief returns the N matrix for IGA elements at a parameter, which fits to the knot vector (e.g. 3D: N & 0 & 0 \\ 0 & N & 0 \\ 0 & 0 & N ...)
    //! @param rCoordinates ... parameter
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const override= 0;

    // -- derivatives --//

    //! @brief returns specific derivative shape functions at an integration point
    //! @param rIP ... id of the integration point
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(int rIP, const Eigen::VectorXi &rKnotIDs) const override= 0;

    //! @brief returns specific derivative shape functions at a parameter, which fits to the knot vector
    //! @param rCoordinates ... parameter
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const override= 0;

    //********************************************
    //       DERIVATIVE SHAPE FUNCTIONS NATURAL
    //********************************************

    //! @brief returns specific derivative shape functions natural via the IP index
    //! @param rIP ... integration point index
    //! @return ... specific derivative shape functions natural
    virtual const Eigen::MatrixXd& GetDerivativeShapeFunctionsNatural(int rIP) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Since the shape functions are calculated on the fly, just use 'GetDerivativeShapeFunctionsNatural' routine!");
    }

    //! @brief returns specific derivative shape functions natural via coordinates
    //! @param rCoordinates ... integration point coordinates
    //! @return ... specific derivative shape functions natural
    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override = 0;

    //********************************************
    //       SURFACE PARAMETRIZATION
    //********************************************

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Use the function 'CalculateNaturalSurfaceCoordinates(rNaturalSurfaceCoordinates, rSurface, rKnots)' instead!");
    }

    //! @brief returns the natural coordinates of the elements surface
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... natural coordinates of the elements surface
    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface, const Eigen::MatrixXd &rKnots) const override = 0;

    //! @brief returns the derivative of the surface parametrization
    //! @param rNaturalSurfaceCoordinates ... natural surface coordinates
    //! @param rSurface ... index of the surface, see documentation of the specific InterpolationType
    //! @return ... derivative of the surface parametrization
    virtual Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override = 0;

    //! @brief returns the number of surfaces
    virtual int GetNumSurfaces() const override = 0;

    //! @brief return the local dimension of the interpolation
    virtual int GetLocalDimension() const override = 0;

    virtual Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const override = 0;

    int GetSurfaceDegree(int rSurface) const override = 0;


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
    virtual std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override = 0;


    //! @brief returns true if a node is on the surface
    //! @param rSurface ... surface id
    //! @param rNaturalNodeCoordinate ... natural coordinate of the node to test
    virtual bool NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "No natural node coordinates in IGA!");
    }

    //! @brief stores the integration point coordinates
    void UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType) override = 0;

    //! @brief return the number node depending the shape and the order
    virtual int CalculateNumNodes() const override = 0;

    //! @brief this method sets the mNumDofs, mNumNodes and mNodeIndices members
    //! @remark it should be called from the ctor InterpolationTypeBase()
    //! but uses pure virutal functions. Thus it must be called in the
    //! ctors of the child classes.
    void Initialize();

    //! @brief transforms unit interval [-1, 1] to the interval [firstKnotCoordinate, secondKnotCoordinate]
    inline double transformation(double rIPCoordinate, double firstKnotCoordinate, double secondKnotCoordinate) const
    {
        return (firstKnotCoordinate + 0.5*(rIPCoordinate + 1)*(secondKnotCoordinate - firstKnotCoordinate));
    }
};
} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::InterpolationBaseIGA)
#endif


