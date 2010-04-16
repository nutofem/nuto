/*******************************************************************************
 Bauhaus-Universitaet Weimar
 Author: Joerg F. Unger,  September 2009
*******************************************************************************/
#include <stdio.h>
#include "nuto/metamodel/TransferFunction.h"
#include "nuto/metamodel/MetamodelException.h"
#include "math.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(NuTo::EmptyTransferFunction,"EmptyTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::HardLimTransferFunction,"HardLimTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::HardLimsTransferFunction,"HardLimsTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::PureLinTransferFunction,"PureLinTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::SatLinTransferFunction,"SatLinTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::SatLinsTransferFunction,"SatLinsTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::LogSigTransferFunction,"LogSigTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::TanSigTransferFunction,"TanSigTransferFunction")
BOOST_CLASS_EXPORT_GUID(NuTo::PosLinTransferFunction,"PosLinTransferFunction")
#endif // ENABLE_SERIALIZATION

namespace NuTo
{
double EmptyTransferFunction::evaluate(double x)
{
    throw MetamodelException("EmptyTransferFunction::evaluate : trying to evaluate empty activation function.");
}

double EmptyTransferFunction::derivative(double x)
{
    throw MetamodelException("EmptyTransferFunction::derivative : trying to evaluate empty activation function.");
}

double EmptyTransferFunction::second_derivative(double x)
{
    throw MetamodelException("EmptyTransferFunction::second_derivative : trying to evaluate empty activation function.");
}

TransferFunction* EmptyTransferFunction::clone()const
{
    return new EmptyTransferFunction();
}

void EmptyTransferFunction::info()const
{
    printf("no activation function\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void EmptyTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void EmptyTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void EmptyTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void EmptyTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void EmptyTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void EmptyTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void EmptyTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double HardLimTransferFunction::evaluate(double x)
{
    if (x<0)
        return 0.;
    else
        return 1.;
}

double HardLimTransferFunction::derivative(double x)
{
    return 0.;
}

double HardLimTransferFunction::second_derivative(double x)
{
    return 0.;
}

TransferFunction* HardLimTransferFunction::clone()const
{
    return new HardLimTransferFunction();
}

void HardLimTransferFunction::info()const
{
    printf("activation function : hardlim (x<0 : 0 else 1)\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void HardLimTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void HardLimTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void HardLimTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void HardLimTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void HardLimTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void HardLimTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void HardLimTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double HardLimsTransferFunction::evaluate(double x)
{
    if (x<0)
        return -1.;
    else
        return 1.;
}

double HardLimsTransferFunction::derivative(double x)
{
    return 0;
}

double HardLimsTransferFunction::second_derivative(double x)
{
    return 0;
}

TransferFunction* HardLimsTransferFunction::clone()const
{
    return new HardLimsTransferFunction();
}

void HardLimsTransferFunction::info()const
{
    printf("activation function : hardlims (x<0 : -1 else 1)\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void HardLimsTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void HardLimsTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void HardLimsTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void HardLimsTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void HardLimsTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void HardLimsTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void HardLimsTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double PureLinTransferFunction::evaluate(double x)
{
    return x;
}

double PureLinTransferFunction::derivative(double x)
{
    return 1;
}

double PureLinTransferFunction::second_derivative(double x)
{
    return 0;
}


TransferFunction* PureLinTransferFunction::clone()const
{
    return new PureLinTransferFunction();
}

void PureLinTransferFunction::info()const
{
    printf("activation function : purelin (x)\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void PureLinTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void PureLinTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void PureLinTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void PureLinTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void PureLinTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void PureLinTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void PureLinTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double SatLinTransferFunction::evaluate(double x)
{
    if (x<0)
        return 0.;
    else if (x<1.)
        return x;
    else
        return 1.;
}

double SatLinTransferFunction::derivative(double x)
{
    if (x<0)
        return 0.;
    else if (x<1.)
        return 1.;
    else
        return 0.;
}

double SatLinTransferFunction::second_derivative(double x)
{
    return 0.;
}

TransferFunction* SatLinTransferFunction::clone()const
{
    return new SatLinTransferFunction();
}

void SatLinTransferFunction::info()const
{
    printf("activation function : satlin (x<0 :0; x<1 x; else 1)\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void SatLinTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void SatLinTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void SatLinTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void SatLinTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void SatLinTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void SatLinTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void SatLinTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double SatLinsTransferFunction::evaluate(double x)
{
    if (x<-1)
        return -1.;
    else if (x<1.)
        return x;
    else
        return 1.;
}

double SatLinsTransferFunction::derivative(double x)
{
    if (x<-1)
        return 0.;
    else if (x<1.)
        return 1;
    else
        return 0.;
}

double SatLinsTransferFunction::second_derivative(double x)
{
    return 0.;
}

TransferFunction* SatLinsTransferFunction::clone()const
{
    return new SatLinsTransferFunction();
}

void SatLinsTransferFunction::info()const
{
    printf("activation function : satlins (x<-1 :-1; x<1 x; else 1)\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void SatLinsTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void SatLinsTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void SatLinsTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void SatLinsTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void SatLinsTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void SatLinsTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void SatLinsTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double LogSigTransferFunction::evaluate(double x)
{
    return (1./(1.+exp(-x)));
}

double LogSigTransferFunction::derivative(double x)
{
    double f = evaluate(x);
	return f*(1.-f);
}

double LogSigTransferFunction::second_derivative(double x)
{
    double f = evaluate(x);
	return f*(1.-f)*(1.-2.*f);
}

TransferFunction* LogSigTransferFunction::clone()const
{
    return new LogSigTransferFunction();
}

void LogSigTransferFunction::info()const
{
    printf("activation function : logsig (1/(1+e^-x))\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void LogSigTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void LogSigTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void LogSigTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void LogSigTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void LogSigTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void LogSigTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void LogSigTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double TanSigTransferFunction::evaluate(double x)
{
    x*=2./3.;
    if (fabs(x)<20)
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

double TanSigTransferFunction::derivative(double x)
{
    x*=2./3.;
    if (fabs(x)<20)
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

double TanSigTransferFunction::second_derivative(double x)
{
    x*=2./3.;
    if (fabs(x)<20)
    {
        double epx = exp(x);
        double emx = exp(-x);
        double value = (epx+emx);
        return -1.7159*32.*(epx-emx)/(9.*value*value*value);
    }
    else
        return 0.;
}

TransferFunction* TanSigTransferFunction::clone()const
{
    return new TanSigTransferFunction();
}
void TanSigTransferFunction::info()const
{
    printf("activation function : tansig ((e^x+e^-x)/(e^x+e^-x))\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void TanSigTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void TanSigTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void TanSigTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void TanSigTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void TanSigTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void TanSigTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void TanSigTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION

double PosLinTransferFunction::evaluate(double x)
{
    if (x<0)
        return 0.;
    else
        return x;
}

double PosLinTransferFunction::derivative(double x)
{
    if (x<0)
        return 0.;
    else
        return 1.;
}

double PosLinTransferFunction::second_derivative(double x)
{
    return 0.;
}

TransferFunction* PosLinTransferFunction::clone()const
{
    return new PosLinTransferFunction();
}

void PosLinTransferFunction::info()const
{
    printf("activation function : poslin (x<0 :0 else x)\n");
    return ;
}

#ifdef ENABLE_SERIALIZATION
template void PosLinTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void PosLinTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void PosLinTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void PosLinTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void PosLinTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void PosLinTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void PosLinTransferFunction::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
}
#endif // ENABLE_SERIALIZATION
} // namespace surrogate

