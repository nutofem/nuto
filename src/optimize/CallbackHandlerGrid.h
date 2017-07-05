#pragma once

#include <boost/dynamic_bitset.hpp>
#include "optimize/OptimizeException.h"
#include <vector>

namespace NuTo
{
//! @author Andrea Keszler, ISM
//! @date July 2010
//! @brief ... abstract class to handle callback routines
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


    virtual std::vector<double>&  GetParameters()
    {
		throw OptimizeException("[CallbackHandlerGrid::SetParameters] SetParameters function not implemented in CallbackHandlerGrid object.");
    }

    virtual std::vector<double>&  GetRightHandSide()
    {
		throw OptimizeException("[CallbackHandlerGrid::GetRightHandSide] GetRightHandSide function not implemented in CallbackHandlerGrid object.");
    }
    virtual void SetParameters(std::vector<double>& rParameters)
    {
		throw OptimizeException("[CallbackHandlerGrid::SetParameters] SetParameters function not implemented in CallbackHandlerGrid object.");
    }

    virtual void SetRightHandSide(std::vector<double>& rRightHandSide)
    {
		throw OptimizeException("[CallbackHandlerGrid::SetRightHandSide] SetRightHandSide function not implemented in CallbackHandlerGrid object.");
    }

    virtual void Gradient (std::vector<double>& rValue,std::vector<double>& rGradient)
	{
		throw OptimizeException("[CallbackHandlerGrid::Gradient] Gradient function not implemented in CallbackHandlerGrid object.");
	}

 	virtual void Hessian (std::vector<double>& rHessian)
	{
		throw OptimizeException("[CallbackHandlerGrid::Gradient] Gradient function not implemented in CallbackHandlerGrid object.");
	}

	//! @brief get DisplacementConstaints
	//! @return dynamic_bitset of constraints
	virtual const boost::dynamic_bitset<> GetDisplacementConstaints()
	{
		throw OptimizeException("[CallbackHandlerGrid::GetDisplacementConstaints] GetDisplacementConstaints function not implemented in CallbackHandlerGrid object.");
	}

	//! @brief correct solution for hanging nodes
	//! @param displacement solution
	virtual void HangingNodesCorrection(std::vector<double>& u)
	{
		throw OptimizeException("[CallbackHandlerGrid::HangingNodesCorrection] HangingNodesCorrection function not implemented in CallbackHandlerGrid object.");
	}


	virtual void SetMisesWielandt (bool rMisesWielandt)
	{
		throw OptimizeException("[CallbackHandlerGrid::SetMisesWielandt] SetMisesWielandt function not implemented in CallbackHandlerGrid object.");
	}
	virtual double GetWeightingFactor()
	{
		throw OptimizeException("[CallbackHandlerGrid::GetWeightingFactor] GetWeightingFactor function not implemented in CallbackHandlerGrid object.");
	}
	virtual void SetWeightingFactor(double rWeight)
	{
		throw OptimizeException("[CallbackHandlerGrid::GetWeightingFactor] GetWeightingFactor function not implemented in CallbackHandlerGrid object.");

	}

private:
    const std::vector<double> *mCallbackSetParameters;
    std::vector<double> *mCallbackGradient;
    std::vector<double> *mCallbackHessian;
};
} //namespace NuTo

