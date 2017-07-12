#pragma once

#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

namespace NuTo
{

class StructureBase;

namespace IpData
{
enum class eIpStaticDataType;
} // namespace IpData

class ResultElementIpData : public ResultBase
{
public:
    //! @brief constructor
    //! @param rFileName:   file name
    //! @param rElementId:  element id
    //! @param rIpDataType: data type at the integration points, e.g. stress, strain, damage...
    ResultElementIpData(const std::string& rFileName, int rElementId, NuTo::IpData::eIpStaticDataType rIpDataType);

    //! @brief calculates the relevant integration point data and adds them to the internal routine
    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot);

    //! @brief calculates the relevant integration point data
    void CalculateValues(const StructureBase& rStructure, Eigen::Matrix<double, 1, Eigen::Dynamic>& rValues) const;

    //! @brief number of data points per time step, e.g. number of stress components for an integration point
    int GetNumData(const StructureBase& rStructure) const override;

    NuTo::eTimeIntegrationResultType GetResultType() const override;

    //! @brief returns the class name
    std::string GetTypeId() const
    {
        return std::string("ResultElementIpValue");
    }

    //! @brief this is used to cast an object of ResultBase to an object of ResultElementIpData
    ResultElementIpData* AsResultElementIpData() override
    {
        return this;
    }

    //! @brief ... Info routine that prints general information about the object
    void Info() const override;

    std::unique_ptr<ResultBase> Clone() const override
    {
        return std::make_unique<ResultElementIpData>(*this);
    }

private:
    int mElementId;
    NuTo::IpData::eIpStaticDataType mIpDataType;
};
} // namespace NuTo
