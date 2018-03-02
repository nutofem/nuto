#pragma once

#include <boost/dynamic_bitset.hpp>
#include "base/Exception.h"
#include <vector>

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief abstract class to handle callback routines
class CallbackHandlerGrid
{

public:
    CallbackHandlerGrid()
    {
        mCallbackSetParameters = nullptr;
        mCallbackGradient = nullptr;
        mCallbackHessian = nullptr;
    }

    virtual ~CallbackHandlerGrid() = default;


    virtual std::vector<double>& GetParameters()
    {
        throw Exception("[CallbackHandlerGrid::SetParameters] SetParameters function not implemented in "
                        "CallbackHandlerGrid object.");
    }

    virtual std::vector<double>& GetRightHandSide()
    {
        throw Exception("[CallbackHandlerGrid::GetRightHandSide] GetRightHandSide function not implemented in "
                        "CallbackHandlerGrid object.");
    }
    virtual void SetParameters(std::vector<double>&)
    {
        throw Exception("[CallbackHandlerGrid::SetParameters] SetParameters function not implemented in "
                        "CallbackHandlerGrid object.");
    }

    virtual void SetRightHandSide(std::vector<double>&)
    {
        throw Exception("[CallbackHandlerGrid::SetRightHandSide] SetRightHandSide function not implemented in "
                        "CallbackHandlerGrid object.");
    }

    virtual void Gradient(std::vector<double>&, std::vector<double>&)
    {
        throw Exception(
                "[CallbackHandlerGrid::Gradient] Gradient function not implemented in CallbackHandlerGrid object.");
    }

    virtual void Hessian(std::vector<double>&)
    {
        throw Exception(
                "[CallbackHandlerGrid::Gradient] Gradient function not implemented in CallbackHandlerGrid object.");
    }

    //! @brief get DisplacementConstaints
    //! @return dynamic_bitset of constraints
    virtual const boost::dynamic_bitset<> GetDisplacementConstaints()
    {
        throw Exception("[CallbackHandlerGrid::GetDisplacementConstaints] GetDisplacementConstaints function "
                        "not implemented in CallbackHandlerGrid object.");
    }

    //! @brief correct solution for hanging nodes
    virtual void HangingNodesCorrection(std::vector<double>&)
    {
        throw Exception("[CallbackHandlerGrid::HangingNodesCorrection] HangingNodesCorrection function not "
                        "implemented in CallbackHandlerGrid object.");
    }


    virtual void SetMisesWielandt(bool)
    {
        throw Exception("[CallbackHandlerGrid::SetMisesWielandt] SetMisesWielandt function not implemented in "
                        "CallbackHandlerGrid object.");
    }
    virtual double GetWeightingFactor()
    {
        throw Exception("[CallbackHandlerGrid::GetWeightingFactor] GetWeightingFactor function not implemented "
                        "in CallbackHandlerGrid object.");
    }
    virtual void SetWeightingFactor(double)
    {
        throw Exception("[CallbackHandlerGrid::GetWeightingFactor] GetWeightingFactor function not implemented "
                        "in CallbackHandlerGrid object.");
    }

private:
    const std::vector<double>* mCallbackSetParameters;
    std::vector<double>* mCallbackGradient;
    std::vector<double>* mCallbackHessian;
};
} // namespace NuTo
