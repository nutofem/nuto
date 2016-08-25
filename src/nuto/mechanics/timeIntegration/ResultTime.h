// $Id: $

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard result class for time
class ResultTime : public ResultBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    ResultTime(const std::string& rIdent);

    std::string GetTypeId() const
    {
    	return std::string("ResultTime");
    }

    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot, double rTime);

    NuTo::eTimeIntegrationResultType GetResultType()const;

    //! @brief number of data points per time step (e.g. number of displacement components of a node)
    int GetNumData(const StructureBase& rStructure)const
    {
    	return 1;
    }

    ResultTime* AsResultTime()
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
};
}

//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultTime)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
