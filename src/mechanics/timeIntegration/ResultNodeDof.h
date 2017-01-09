// $Id: $

#pragma once


#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/timeIntegration/ResultBase.h"


namespace NuTo
{

class StructureBase;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
class NodeBase;
class ResultNodeDof : public ResultBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    ResultNodeDof(const std::string& rIdent, int rNodeId);

    //! @brief calculate the relevant nodal dofs and add to the internal routine
    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot);

    //! @brief calculate the relevant nodal dofs
    virtual void CalculateValues(const StructureBase& rStructure, Eigen::Matrix<double, 1, Eigen::Dynamic>& rValues)const=0;

    ResultNodeDof* AsResultNodeDof()
    {
    	return this;
    }

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
    int mNodeId;
};
}

//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultNodeDof)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
