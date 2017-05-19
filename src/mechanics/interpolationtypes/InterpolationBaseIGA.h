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

    //********************************************
    //             NODE METHODS
    //********************************************

    const Eigen::VectorXd& GetNaturalNodeCoordinates(int rNodeIndex) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "No natural node coordinates in IGA!");
    }

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

    virtual const Eigen::VectorXd& ShapeFunctions(const Eigen::VectorXd& naturalCoordinates) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "The shape functions are calculated on the fly, use 'GetShapeFunctions' routine!");
    }

    virtual const Eigen::MatrixXd& MatrixN(const Eigen::VectorXd& naturalCoordinates) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "The shape functions are calculated on the fly, use 'GetMatrixN' routine!");
    }

    //********************************************
    //       DERIVATIVE SHAPE FUNCTIONS NATURAL
    //********************************************

    virtual const Eigen::MatrixXd& DerivativeShapeFunctionsNatural(const Eigen::VectorXd& naturalCoordinates) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Since the shape functions are calculated on the fly, just use 'GetDerivativeShapeFunctionsNatural' routine!");
    }

    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override = 0;

    //********************************************
    //       SURFACE PARAMETRIZATION
    //********************************************

    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Use the function 'CalculateNaturalSurfaceCoordinates(rNaturalSurfaceCoordinates, rSurface, rKnots)' instead!");
    }

    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface, const Eigen::MatrixXd &rKnots) const override = 0;


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

    virtual std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override = 0;

    virtual bool NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "No natural node coordinates in IGA!");
    }

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


