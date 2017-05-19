/*
 * InterpolationBaseFEM.h
 *
 *  Created on: 17 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/MechanicsException.h"

#include <vector>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "math/CustomBoostSerializationExtensions.h"
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

    const Eigen::VectorXd& GetNaturalNodeCoordinates(int rNodeIndex) const override;

    virtual Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndex) const override = 0;

    void CalculateSurfaceNodeIds() override;

    //********************************************
    //       SHAPE FUNCTIONS
    //********************************************

    const Eigen::VectorXd& GetShapeFunctions(int rIP) const override;

    const Eigen::MatrixXd& GetMatrixN(int rIP) const override;

    virtual Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override = 0;

    Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const override;


    // --- IGA interpolation--- //

    virtual Eigen::VectorXd CalculateShapeFunctions(int rIP, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    virtual Eigen::MatrixXd CalculateMatrixN(int rIP, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    virtual Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    virtual Eigen::MatrixXd CalculateMatrixNDerivative(const Eigen::VectorXd& rParameters, const Eigen::VectorXi& rKnotIDs, int rDerivative, int rDirection) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "So far implemeneted only for IGA!");
    }
    //********************************************
    //       DERIVATIVE SHAPE FUNCTIONS NATURAL
    //********************************************

    const Eigen::MatrixXd & GetDerivativeShapeFunctionsNatural(int rIP) const override;

    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override = 0;

    // --- IGA interpolation--- //

    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(int rIP, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    virtual Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    //********************************************
    //       SURFACE PARAMETRIZATION
    //********************************************

    virtual Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override = 0;

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface, const Eigen::MatrixXd &rKnots) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    virtual Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override = 0;

    virtual int GetNumSurfaces() const override = 0;

    virtual int GetLocalDimension() const override = 0;

    Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "IGA specific function!");
    }

    int GetSurfaceDegree(int rSurface) const override
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

    //! @brief calculate and store the shape functions and their derivatives
    //! @param rIntegrationType ... integration type
    void UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType) override;

protected:

    virtual std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override = 0;

    bool NodeIsOnSurface(int rSurface, const Eigen::VectorXd& rNaturalNodeCoordinate) const override;

    virtual int CalculateNumNodes() const override = 0;

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
