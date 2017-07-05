// $Id: $
#pragma once

#include "mechanics/timeIntegration/ResultBase.h"

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
public:
    //! @brief constructor
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

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

protected:
    int mGroupNodeId;
};
}

// namespace NuTo
