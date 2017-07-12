#pragma once

#include <eigen3/Eigen/Core>
#include <memory>
#include "mechanics/MechanicsException.h"

namespace NuTo
{
class ResultNodeDof;
class ResultGroupNodeDof;
class ResultTime;
class ResultElementIpData;
class StructureBase;
enum class eTimeIntegrationResultType;

//! @author Jörg F. Unger, BAM
//! @date December 2013
//! @brief ... standard abstract class for all results
class ResultBase
{
public:
    //! @brief constructor
    ResultBase(const std::string& rIdent);

    //! @brief deconstructor
    virtual ~ResultBase();

    virtual std::unique_ptr<ResultBase> Clone() const = 0;

    void SetIdent(const std::string& rIdent);

    std::string GetIdent() const;

    void Resize(const StructureBase& rStructure, int rNumResultSteps, bool rInitialize);

    void WriteToFile(const std::string& rResultDir, int rNumTimeSteps) const;

    //! @brief number of data points per time step (e.g. number of displacement components of a node)
    virtual int GetNumData(const StructureBase& rStructure) const = 0;

    virtual NuTo::eTimeIntegrationResultType GetResultType() const = 0;

    virtual ResultNodeDof* AsResultNodeDof()
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "object is not of this type.");
    }

    virtual ResultTime* AsResultTime()
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "object is not of this type.");
    }

    virtual ResultGroupNodeDof* AsResultGroupNodeDof()
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "object is not of this type.");
    }


    virtual ResultElementIpData* AsResultElementIpData()
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "object is not of this type.");
    }

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    virtual void Info() const
    {
    }

protected:
    std::string mIdent;
    Eigen::MatrixXd mData;
};

ResultBase* new_clone(const ResultBase& result);
} // namespace NuTo
