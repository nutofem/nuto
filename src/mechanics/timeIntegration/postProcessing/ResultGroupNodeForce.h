#pragma once

#include "mechanics/timeIntegration/postProcessing/ResultGroupNodeDof.h"

namespace NuTo
{

class ResultGroupNodeForce : public ResultGroupNodeDof
{
public:
    ResultGroupNodeForce(const std::string& rIdent, int rGroupNodeId);

    //! @brief number of dofs (e.g. number of displacement components of a node
    int GetNumData(const StructureBase& rStructure) const override;

    NuTo::eTimeIntegrationResultType GetResultType() const override;

    Eigen::VectorXd CalculateValues(const StructureBase& rStructure, const Eigen::VectorXd& rResidual_j,
                                    const Eigen::VectorXd& rResidual_k) const override;

    std::string GetTypeId() const
    {
        return std::string("ResultGroupNodeForce");
    }

    void Info() const override
    {
    }

    std::unique_ptr<ResultBase> Clone() const override
    {
        return std::make_unique<ResultGroupNodeForce>(*this);
    }
};
} // namespace NuTo
