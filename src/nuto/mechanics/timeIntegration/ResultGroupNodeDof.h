// $Id: $

#ifndef ResultGroupNodeDof_H
#define ResultGroupNodeDof_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultBase.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
template <class T>
class Group;

class NodeBase;

class ResultGroupNodeDof : public ResultBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    ResultGroupNodeDof(const std::string& rIdent, int rNodeGroupId);

    NuTo::ResultGroupNodeDof* AsResultGroupNodeDof()override
    {
    	return this;
    }

    virtual void CalculateValues(const StructureBase& rStructure,
    		   const FullVector<double,Eigen::Dynamic>& rResidual_j, const FullVector<double,Eigen::Dynamic>& rResidual_k,
    		   FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rResult)const=0;

    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot,
    		const FullVector<double,Eigen::Dynamic>& rResidual_j,
    		const FullVector<double,Eigen::Dynamic>& rResidual_k);

    const NuTo::Group<NodeBase>* GetGroupNodePtr(const StructureBase& rStructure)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

protected:
    int mGroupNodeId;
};
}

//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultGroupNodeDof)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // ResultGroupNodeDof_H
