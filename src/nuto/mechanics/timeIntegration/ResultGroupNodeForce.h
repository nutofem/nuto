// $Id: $

#ifndef ResultGroupNodeForce_H
#define ResultGroupNodeForce_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultGroupNodeDof.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

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
    int GetNumData(const StructureBase& rStructure)const;

    NuTo::TimeIntegration::eResultType GetResultType()const
    {
    	return NuTo::TimeIntegration::GROUP_NODE_FORCE;
    }

    void CalculateValues(const StructureBase& rStructure,
    		const FullVector<double,Eigen::Dynamic>& rResidual_j,
    		const FullVector<double,Eigen::Dynamic>& rResidual_k,
    		FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)const;

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
    void Info()const
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
#endif // ResultGroupNodeForce_H
