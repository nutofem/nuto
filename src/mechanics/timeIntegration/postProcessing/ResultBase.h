#pragma once

#include <eigen3/Eigen/Core>
#include <memory>
#include "mechanics/MechanicsException.h"

namespace NuTo
{
class StructureBase;
class StructureOutputBlockVector;

//! @brief Abstract class for all results
class ResultBase
{
public:
    //! @brief constructor
    ResultBase(const std::string& rIdent);

    //! @brief deconstructor
    virtual ~ResultBase();

    virtual std::unique_ptr<ResultBase> Clone() const = 0;

    void Resize(const StructureBase& rStructure, int rNumResultSteps, bool rInitialize);

    virtual void CalculateAndAddValues(const StructureBase& rStructure, int timeStep,
                                       const StructureOutputBlockVector& residual, double currentTime) = 0;

    void WriteToFile(const std::string& rResultDir, int rNumTimeSteps) const;

    //! @brief number of data points per time step (e.g. number of displacement components of a node)
    virtual int GetNumData(const StructureBase& rStructure) const = 0;

protected:
    std::string mIdent;
    Eigen::MatrixXd mData;
};

ResultBase* new_clone(const ResultBase& result);
} // namespace NuTo
