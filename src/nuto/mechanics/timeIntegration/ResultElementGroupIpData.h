#pragma once


#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultBase.h"


namespace NuTo
{

class StructureBase;

template <class T, int rows, int cols> class  FullMatrix;

namespace IpData
{
    enum class eIpStaticDataType;
}// namespace IpData

//! @author Peter Otto
//! @date Dec 2016
//! @brief Outputs integration point values with its coordinates for a whole element group
class ResultElementGroupIpData : public ResultBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    //! @param rFileName:   file name
    //! @param rElementId:  element id
    //! @param rIpDataType: data type at the integration points, e.g. stress, strain, damage...
    ResultElementGroupIpData(const std::string& rIdent, int rElementGroupId, int rComponent, const std::vector<int> &rIP, NuTo::IpData::eIpStaticDataType rIpDataType);

    //! @brief calculates the relevant integration point data and adds them to the internal routine
    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot);

    //! @brief calculates the relevant integration point data
    void CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rValues)const;

    int GetNumData(const StructureBase& rStructure)const
    {
        throw MechanicsException(std::string(__PRETTY_FUNCTION__) + "\t: Use 'GetNumData(const StructureBase& rStructure, int &rows, int &cols)' instead!.");
    }

    //! @brief number of data points per time step, e.g. number of stress components for an integration point
    void GetNumData(const StructureBase& rStructure, int &rows, int &cols) const override;

    NuTo::eTimeIntegrationResultType GetResultType() const;

    //! @brief returns the class name
    std::string GetTypeId() const
    {
    	return std::string("ResultElementIpValue");
    }

    //! @brief this is used to cast an object of ResultBase to an object of ResultElementGroupIpData
    ResultElementGroupIpData* AsResultElementGroupIpData()
    {
        return this;
    }

    //! @brief ... Info routine that prints general information about the object
    void Info()const;

private:
    //! @brief the group id
    int mElementGroupId;
    //! @brief the component of the static data type
    int mComponent;
    //! @brief the integration point ids to write out
    std::vector<int> mIP;
    //! @brief the static data type (stress, straint, etc.)
    NuTo::IpData::eIpStaticDataType mIpDataType;
};
}// namespace NuTo

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ResultElementIpValue)
#endif      // SWIG
#endif      // ENABLE_SERIALIZATION
