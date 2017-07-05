// $Id: $

#pragma once

#include "mechanics/timeIntegration/ResultNodeDof.h"


namespace NuTo
{

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all results
class ResultNodeAcceleration : public ResultNodeDof
{
public:
    //! @brief constructor
    ResultNodeAcceleration(const std::string& rIdent, int rNodeId);

    //! @brief calculate the relevant nodal dofs
    Eigen::VectorXd CalculateValues(const StructureBase& rStructure) const override;

    //! @brief number of data points per time step (e.g. number of Accelerationlacement components of a node
    int GetNumData(const StructureBase& rStructure)const override;

    NuTo::eTimeIntegrationResultType GetResultType()const override;

    std::string GetTypeId() const
    {
    	return std::string("ResultNodeAcceleration");
    }

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override
    {

    }

protected:
};
}

//namespace NuTo
