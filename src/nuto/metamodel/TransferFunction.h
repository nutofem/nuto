/*******************************************************************************
 Bauhaus-Universitaet Weimar
 Author: Joerg F. Unger,  September 2009
*******************************************************************************/

#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/assume_abstract.hpp>
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
    void serialize(Archive & ar, const unsigned int version){}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const {return EMPTY_TRANSFER_FUNCTION;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const {return HardLim;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const {return HardLimS;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const {return PURELIN;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION
	TransferFunction::eTransferFunction get_enum()const {return SATLIN;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
	TransferFunction::eTransferFunction get_enum()const {return SATLINS;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const {return LOGSIG;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const {return TANSIG;}
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

    double evaluate(double x);
    TransferFunction* clone()const;
    double derivative(double x);
    double second_derivative(double x);
    void info()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class 
    //! @param ar         archive
    //! @param version    version
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

	TransferFunction::eTransferFunction get_enum()const {return POSLIN;}
};

}  // namespace surrogate
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::TransferFunction)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
#endif // TRANSFER_FUNCTION_H
