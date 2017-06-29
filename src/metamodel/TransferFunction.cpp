/*******************************************************************************
 Bauhaus-Universitaet Weimar
 Author: Joerg F. Unger,  September 2009
*******************************************************************************/
#include "metamodel/TransferFunction.h"
#include "base/Exception.h"
#include <cmath>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
//BOOST_CLASS_EXPORT_GUID(NuTo::EmptyTransferFunction,"EmptyTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::HardLimTransferFunction,"HardLimTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::HardLimsTransferFunction,"HardLimsTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::PureLinTransferFunction,"PureLinTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::SatLinTransferFunction,"SatLinTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::SatLinsTransferFunction,"SatLinsTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::LogSigTransferFunction,"LogSigTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::TanSigTransferFunction,"TanSigTransferFunction")
//BOOST_CLASS_EXPORT_GUID(NuTo::PosLinTransferFunction,"PosLinTransferFunction")
#endif // ENABLE_SERIALIZATION

#ifdef ENABLE_SERIALIZATION
template void NuTo::TransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::TransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::TransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::TransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::TransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::TransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::TransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize TransferFunction" << std::endl;
#endif

#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize TransferFunction" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION

double NuTo::EmptyTransferFunction::evaluate(double x)
{
    throw Exception("EmptyTransferFunction::evaluate : trying to evaluate empty activation function.");
}

double NuTo::EmptyTransferFunction::derivative(double x)
{
    throw Exception("EmptyTransferFunction::derivative : trying to evaluate empty activation function.");
}

double NuTo::EmptyTransferFunction::second_derivative(double x)
{
    throw Exception("EmptyTransferFunction::second_derivative : trying to evaluate empty activation function.");
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

#ifdef ENABLE_SERIALIZATION
template void NuTo::EmptyTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EmptyTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EmptyTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EmptyTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EmptyTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::EmptyTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::EmptyTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EmptyTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EmptyTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::EmptyTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::HardLimTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::HardLimTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::HardLimTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::HardLimTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::HardLimTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::HardLimTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::HardLimTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize HardLimTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize HardLimTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::HardLimTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::HardLimsTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::HardLimsTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::HardLimsTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::HardLimsTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::HardLimsTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::HardLimsTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::HardLimsTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize HardLimsTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize HardLimsTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::HardLimsTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::PureLinTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::PureLinTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::PureLinTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::PureLinTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::PureLinTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::PureLinTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::PureLinTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize PureLinTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize PureLinTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::PureLinTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::SatLinTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SatLinTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SatLinTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SatLinTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SatLinTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::SatLinTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::SatLinTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SatLinTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SatLinTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SatLinTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::SatLinsTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SatLinsTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SatLinsTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SatLinsTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SatLinsTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::SatLinsTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::SatLinsTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SatLinsTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SatLinsTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SatLinsTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::LogSigTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::LogSigTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::LogSigTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::LogSigTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::LogSigTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::LogSigTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::LogSigTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LogSigTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LogSigTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LogSigTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::TanSigTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::TanSigTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::TanSigTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::TanSigTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::TanSigTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::TanSigTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::TanSigTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize TanSigTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize TanSigTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::TanSigTransferFunction)
#endif // ENABLE_SERIALIZATION

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

#ifdef ENABLE_SERIALIZATION
template void NuTo::PosLinTransferFunction::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::PosLinTransferFunction::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::PosLinTransferFunction::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::PosLinTransferFunction::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::PosLinTransferFunction::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::PosLinTransferFunction::serialize(boost::archive::text_oarchive & ar, const unsigned int version);

template<class Archive>
void NuTo::PosLinTransferFunction::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize PosLinTransferFunction" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TransferFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize PosLinTransferFunction" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::PosLinTransferFunction)
#endif // ENABLE_SERIALIZATION

