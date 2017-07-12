#pragma once

#include "mechanics/timeIntegration/postProcessing/ResultBase.h"

namespace NuTo
{
template <class T>
class Group;

class NodeBase;

class ResultGroupNodeDof : public ResultBase
{
public:
    ResultGroupNodeDof(const std::string& rIdent, int rNodeGroupId);

    NuTo::ResultGroupNodeDof* AsResultGroupNodeDof() override
    {
        return this;
    }

    virtual Eigen::VectorXd CalculateValues(const StructureBase& rStructure, const Eigen::VectorXd& rResidual_j,
                                            const Eigen::VectorXd& rResidual_k) const = 0;

    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot, const Eigen::VectorXd& rResidual_j,
                               const Eigen::VectorXd& rResidual_k);

    const NuTo::Group<NodeBase>* GetGroupNodePtr(const StructureBase& rStructure) const;

    void Info() const override;

protected:
    int mGroupNodeId;
};
} // namespace NuTo
