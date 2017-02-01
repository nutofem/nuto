// $Id$

/*******************************************************************************
 Bauhaus-Universitaet Weimar
 Author: Joerg F. Unger,  September 2009
*******************************************************************************/

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

namespace NuTo
{

class TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    typedef enum
    {
        EMPTY_TRANSFER_FUNCTION, HardLim, HardLimS, PURELIN, SATLIN, SATLINS,
        LOGSIG, TANSIG, POSLIN
    } eTransferFunction;

    TransferFunction()
    {}

    virtual ~TransferFunction()
    {}

    virtual double evaluate(double x)=0;
    virtual TransferFunction* clone()const=0;
    virtual double derivative(double x)=0;
    virtual double second_derivative(double x)=0;
    virtual void info()const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    virtual eTransferFunction get_enum()const=0;
};

/*******************************************************************************
 none
 *******************************************************************************/
class EmptyTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    EmptyTransferFunction()
    {}

    ~EmptyTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const override {return EMPTY_TRANSFER_FUNCTION;}
};

/*******************************************************************************
 hard limit function
 *******************************************************************************/
class HardLimTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    HardLimTransferFunction()
    {}

    ~HardLimTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const override {return HardLim;}
};

/*******************************************************************************
 symmetric hard limit function
 *******************************************************************************/
class HardLimsTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    HardLimsTransferFunction()
    {}

    ~HardLimsTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const override {return HardLimS;}
};

/*******************************************************************************
 linear limit function
 *******************************************************************************/
class PureLinTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    PureLinTransferFunction()
    {}

    ~PureLinTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const override {return PURELIN;}
};

/*******************************************************************************
 saturating linear limit function
 *******************************************************************************/
class SatLinTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    SatLinTransferFunction()
    {}

    ~SatLinTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION
	TransferFunction::eTransferFunction get_enum()const override {return SATLIN;}
};

/*******************************************************************************
 symmetric saturating linear limit function
 *******************************************************************************/
class SatLinsTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    SatLinsTransferFunction()
    {}

    ~SatLinsTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
	TransferFunction::eTransferFunction get_enum()const override {return SATLINS;}
};

/*******************************************************************************
 log sigmoid function
 *******************************************************************************/
class LogSigTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    LogSigTransferFunction()
    {}

    ~LogSigTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const override {return LOGSIG;}
};

/*******************************************************************************
 hyperbolic tangent sigmoid
 *******************************************************************************/
class TanSigTransferFunction : public TransferFunction
{
public:
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

    TanSigTransferFunction()
    {}

    ~TanSigTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const override {return TANSIG;}
};

/*******************************************************************************
 positive linear
 *******************************************************************************/
class PosLinTransferFunction : public TransferFunction
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    PosLinTransferFunction()
    {}

    ~PosLinTransferFunction()
    {}

    double evaluate(double x) override;
    TransferFunction* clone()const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info()const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const override {return POSLIN;}
};

}  // namespace surrogate
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::TransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::EmptyTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::HardLimTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::HardLimsTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::PureLinTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::SatLinTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::SatLinsTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::LogSigTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::TanSigTransferFunction)
BOOST_CLASS_EXPORT_KEY(NuTo::PosLinTransferFunction)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
