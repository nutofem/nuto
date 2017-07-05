/*******************************************************************************
 Bauhaus-Universitaet Weimar
 Author: Joerg F. Unger,  September 2009
*******************************************************************************/
#include "metamodel/TransferFunction.h"
#include "metamodel/MetamodelException.h"
#include <cmath>

double NuTo::EmptyTransferFunction::evaluate(double x)
{
    throw MetamodelException("EmptyTransferFunction::evaluate : trying to evaluate empty activation function.");
}

double NuTo::EmptyTransferFunction::derivative(double x)
{
    throw MetamodelException("EmptyTransferFunction::derivative : trying to evaluate empty activation function.");
}

double NuTo::EmptyTransferFunction::second_derivative(double x)
{
    throw MetamodelException("EmptyTransferFunction::second_derivative : trying to evaluate empty activation function.");
}

NuTo::TransferFunction* NuTo::EmptyTransferFunction::clone()const
{
    return new EmptyTransferFunction();
}

void NuTo::EmptyTransferFunction::info()const
{
    printf("no activation function\n");
    return ;
}

double NuTo::HardLimTransferFunction::evaluate(double x)
{
    if (x<0)
        return 0.;
    else
        return 1.;
}

double NuTo::HardLimTransferFunction::derivative(double x)
{
    return 0.;
}

double NuTo::HardLimTransferFunction::second_derivative(double x)
{
    return 0.;
}

NuTo::TransferFunction* NuTo::HardLimTransferFunction::clone()const
{
    return new HardLimTransferFunction();
}

void NuTo::HardLimTransferFunction::info()const
{
    printf("activation function : hardlim (x<0 : 0 else 1)\n");
    return ;
}

double NuTo::HardLimsTransferFunction::evaluate(double x)
{
    if (x<0)
        return -1.;
    else
        return 1.;
}

double NuTo::HardLimsTransferFunction::derivative(double x)
{
    return 0;
}

double NuTo::HardLimsTransferFunction::second_derivative(double x)
{
    return 0;
}

NuTo::TransferFunction* NuTo::HardLimsTransferFunction::clone()const
{
    return new HardLimsTransferFunction();
}

void NuTo::HardLimsTransferFunction::info()const
{
    printf("activation function : hardlims (x<0 : -1 else 1)\n");
    return ;
}

double NuTo::PureLinTransferFunction::evaluate(double x)
{
    return x;
}

double NuTo::PureLinTransferFunction::derivative(double x)
{
    return 1;
}

double NuTo::PureLinTransferFunction::second_derivative(double x)
{
    return 0;
}


NuTo::TransferFunction* NuTo::PureLinTransferFunction::clone()const
{
    return new PureLinTransferFunction();
}

void NuTo::PureLinTransferFunction::info()const
{
    printf("activation function : purelin (x)\n");
    return ;
}

double NuTo::SatLinTransferFunction::evaluate(double x)
{
    if (x<0)
        return 0.;
    else if (x<1.)
        return x;
    else
        return 1.;
}

double NuTo::SatLinTransferFunction::derivative(double x)
{
    if (x<0)
        return 0.;
    else if (x<1.)
        return 1.;
    else
        return 0.;
}

double NuTo::SatLinTransferFunction::second_derivative(double x)
{
    return 0.;
}

NuTo::TransferFunction* NuTo::SatLinTransferFunction::clone()const
{
    return new SatLinTransferFunction();
}

void NuTo::SatLinTransferFunction::info()const
{
    printf("activation function : satlin (x<0 :0; x<1 x; else 1)\n");
    return ;
}

double NuTo::SatLinsTransferFunction::evaluate(double x)
{
    if (x<-1)
        return -1.;
    else if (x<1.)
        return x;
    else
        return 1.;
}

double NuTo::SatLinsTransferFunction::derivative(double x)
{
    if (x<-1)
        return 0.;
    else if (x<1.)
        return 1;
    else
        return 0.;
}

double NuTo::SatLinsTransferFunction::second_derivative(double x)
{
    return 0.;
}

NuTo::TransferFunction* NuTo::SatLinsTransferFunction::clone()const
{
    return new SatLinsTransferFunction();
}

void NuTo::SatLinsTransferFunction::info()const
{
    printf("activation function : satlins (x<-1 :-1; x<1 x; else 1)\n");
    return ;
}

double NuTo::LogSigTransferFunction::evaluate(double x)
{
    return (1./(1.+exp(-x)));
}

double NuTo::LogSigTransferFunction::derivative(double x)
{
    double f = evaluate(x);
	return f*(1.-f);
}

double NuTo::LogSigTransferFunction::second_derivative(double x)
{
    double f = evaluate(x);
	return f*(1.-f)*(1.-2.*f);
}

NuTo::TransferFunction* NuTo::LogSigTransferFunction::clone()const
{
    return new LogSigTransferFunction();
}

void NuTo::LogSigTransferFunction::info()const
{
    printf("activation function : logsig (1/(1+e^-x))\n");
    return ;
}

double NuTo::TanSigTransferFunction::evaluate(double x)
{
    x*=2./3.;
    if (std::abs(x)<20)
    {
        double epx = exp(x);
        double emx = exp(-x);
        return 1.7159*(epx-emx)/(epx+emx);
    }
    else
    {
		if (x>0)
            return 1.7159;
        else
            return -1.7159;
    }

}

double NuTo::TanSigTransferFunction::derivative(double x)
{
    x*=2./3.;
    if (std::abs(x)<20)
    {
		double epx = exp(x);
        double emx = exp(-x);
        double value = (epx+emx);
        return 1.7159*8./(3.*value*value);
    }
    else
	{
        return 0.;
	}
}

double NuTo::TanSigTransferFunction::second_derivative(double x)
{
    x*=2./3.;
    if (std::abs(x)<20)
    {
        double epx = exp(x);
        double emx = exp(-x);
        double value = (epx+emx);
        return -1.7159*32.*(epx-emx)/(9.*value*value*value);
    }
    else
        return 0.;
}

NuTo::TransferFunction* NuTo::TanSigTransferFunction::clone()const
{
    return new TanSigTransferFunction();
}
void NuTo::TanSigTransferFunction::info()const
{
    printf("activation function : tansig ((e^x+e^-x)/(e^x+e^-x))\n");
    return ;
}

double NuTo::PosLinTransferFunction::evaluate(double x)
{
    if (x<0)
        return 0.;
    else
        return x;
}

double NuTo::PosLinTransferFunction::derivative(double x)
{
    if (x<0)
        return 0.;
    else
        return 1.;
}

double NuTo::PosLinTransferFunction::second_derivative(double x)
{
    return 0.;
}

NuTo::TransferFunction* NuTo::PosLinTransferFunction::clone()const
{
    return new PosLinTransferFunction();
}

void NuTo::PosLinTransferFunction::info()const
{
    printf("activation function : poslin (x<0 :0 else x)\n");
    return ;
}

