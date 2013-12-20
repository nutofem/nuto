// $Id: $

#ifndef ResultTime_H
#define ResultTime_H

#include <ctime>
#include <array>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultBase.h"
#include "nuto/math/FullVector.h"

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

    void WriteToFile(const std::string& rResultDir, int rTimeStepPlot)const override;

    //! @brief resize the data (to avoid reallocation every time step, a certain number of
    //! @param rNumTimeSteps ... number of time steps to allocate the matrix for
    //! @param rInitValues ... true set all the values to zero, false - leave the values as is
    void Resize(int rNumTimeSteps, bool rInitValues);

    std::string GetTypeId() const
    {
    	return std::string("ResultTime");
    }

    void CalculateAndAddValues(int rTimeStepPlot, double rTime);

    NuTo::TimeIntegration::eResultType GetResultType()const
    {
    	return NuTo::TimeIntegration::TIME;
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
    FullVector<double, Eigen::Dynamic> mData;
};
}

//namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultTime)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION
#endif // ResultTime_H
