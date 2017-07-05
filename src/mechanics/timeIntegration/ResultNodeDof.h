// $Id: $

#pragma once


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
public:
    //! @brief constructor
    ResultNodeDof(const std::string& rIdent, int rNodeId);

    //! @brief calculate the relevant nodal dofs and add to the internal routine
    void CalculateAndAddValues(const StructureBase& rStructure, int rTimeStepPlot);

    //! @brief calculate the relevant nodal dofs
    virtual Eigen::VectorXd CalculateValues(const StructureBase& rStructure) const = 0;

    ResultNodeDof* AsResultNodeDof() override
    {
        return this;
    }

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

protected:
    int mNodeId;
};
}

// namespace NuTo
