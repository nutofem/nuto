// $Id$
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/OptimizeException.h"

NuTo::CallbackHandlerGrid::CallbackHandlerGrid() : CallbackHandler()
{
	mCallbackSetParameters=0;
	mCallbackObjective=0;
	mCallbackGradient=0;
	mCallbackHessian=0;
}
void NuTo::CallbackHandlerGrid::SetParameters(const NuTo::FullMatrix<double>& rParameters)
{
	throw OptimizeException("[CallbackHandler::Objective] Objective function not implemented in CallbackHandler object.");
}

void NuTo::CallbackHandlerGrid::SetParameters(NuTo::FullMatrix<double>& rParameters)
{
	mCallbackSetParameters = &rParameters;
}

double NuTo::CallbackHandlerGrid::Objective()const
{
	throw OptimizeException("[CallbackHandler::Objective] Objective function not implemented in CallbackHandler object.");
	return 0;
}

void NuTo::CallbackHandlerGrid::Gradient (NuTo::FullMatrix<double>& rGradient)const
{
	throw OptimizeException("[CallbackHandler::Gradient] Gradient function not implemented in CallbackHandler object.");
}

void NuTo::CallbackHandlerGrid::Hessian (NuTo::FullMatrix<double>& rHessian)const
{
	throw OptimizeException("[CallbackHandler::Gradient] Gradient function not implemented in CallbackHandler object.");
}

void NuTo::CallbackHandlerGrid::Info()const
{
	std::cout << "CallbackHandlerGrid" << std::endl;
}
