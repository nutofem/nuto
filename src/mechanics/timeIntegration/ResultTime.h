// $Id: $

#pragma once

#include "mechanics/timeIntegration/ResultBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard result class for time
class ResultTime : public ResultBase
{
public:
    //! @brief constructor
    ResultTime(const std::string& rIdent);

    std::string GetTypeId() const
    {
        return std::string("ResultTime");
    }

    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot, double rTime);

    NuTo::eTimeIntegrationResultType GetResultType() const override;

    //! @brief number of data points per time step (e.g. number of displacement components of a node)
    int GetNumData(const StructureBase&) const override
    {
        return 1;
    }

    ResultTime* AsResultTime() override
    {
        return this;
    }

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

protected:
};
}

// namespace NuTo
