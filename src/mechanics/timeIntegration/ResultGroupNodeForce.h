// $Id: $

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/timeIntegration/ResultGroupNodeDof.h"

namespace NuTo
{

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
class ResultGroupNodeForce : public ResultGroupNodeDof
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    ResultGroupNodeForce(const std::string& rIdent, int rGroupNodeId);

    //! @brief number of dofs (e.g. number of displacement components of a node
    int GetNumData(const StructureBase& rStructure) const override;

    NuTo::eTimeIntegrationResultType GetResultType() const override;

    Eigen::VectorXd CalculateValues(const StructureBase& rStructure,
    		const Eigen::VectorXd& rResidual_j,
    		const Eigen::VectorXd& rResidual_k) const override;

    std::string GetTypeId() const
    {
    	return std::string("ResultGroupNodeForce");
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override
    {

    }

protected:
};
}

//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultGroupNodeForce)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
