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

    virtual Eigen::VectorXd CalculateValues(const StructureBase& rStructure,
                                            const StructureOutputBlockVector& residual) const = 0;

    void CalculateAndAddValues(const StructureBase& rStructure, int timeStep,
                               const StructureOutputBlockVector& residual, double currentTime) override;

    const NuTo::Group<NodeBase>* GetGroupNodePtr(const StructureBase& rStructure) const;

    void Info() const override;

protected:
    int mGroupNodeId;
};
} // namespace NuTo
