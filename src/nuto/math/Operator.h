// $Id$

#ifndef OPERATOR_H
#define OPERATOR_H

#include "nuto/math/MathException.h"

namespace NuTo
{

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all Operators in NuTo
template <class T>
class MonadicOperator
{
public:
    MonadicOperator<T>()
    {}
    virtual ~MonadicOperator<T>()
    {}

    //! @brief ... operater with one inputs and one output
    //! @param rInput ... input
    //! @return ... result
    virtual T Evaluate(const T& rInput)const=0;
};
}

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all Operators in NuTo
template <class T>
class DyadicOperator
{
public:
    DyadicOperator<T>()
    {}

    virtual ~DyadicOperator<T>()
    {}

    //! @brief ... operater with two inputs and one output
    //! @param rInput1 ... first input
    //! @param rInput1 ... second input
    //! @return ... result
    virtual T Evaluate(const T& rInput1, const T& rInput2)const=0;
};

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... calculates the absolute value
template <class T>
class MOperatorAbs : public MonadicOperator<T>
{
public:
    MOperatorAbs<T>() :  MonadicOperator<T>()
    {}

    //! @brief ... operater with one input and one output
    //! @param rInput ... input
    //! @return result
    T Evaluate(const T& rInput)const
    {
        if (rInput<0)
            return -rInput;
        return rInput;
    }
};

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... calculates the absolute value
template <class T>
class MOperatorMin : public MonadicOperator<T>
{
public:
    MOperatorMin<T>(const T& rMin) :  MonadicOperator<T>()
    {
        mMin = rMin;
    }

    //! @brief ... operater with one input and one output
    //! @param rInput ... input
    //! @return ... result
    T Evaluate(const T& rInput)const
    {
        if (rInput>mMin)
            return mMin;
        return rInput;
    }

private:
    T mMin;

};
//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... calculates the absolute value
template <class T>
class MOperatorMax : public MonadicOperator<T>
{
public:
    MOperatorMax<T>(const T& rMax) :  MonadicOperator<T>()
    {
        mMax = rMax;
    }

    //! @brief ... operater with one input and one output
    //! @param rInput ... input
    //! @return ... result
    T Evaluate(const T& rInput)const
    {
        if (rInput<mMax)
            return mMax;
        return rInput;
    }

private:
    T mMax;

};

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all Operators in NuTo
template <class T>
class DOperatorMin : public DyadicOperator<T>
{
public:
    DOperatorMin<T>() :  DyadicOperator<T>()
    {}

    //! @brief ... operator with two inputs and one output
    //! @param rInput1 ... first input
    //! @param rInput2 ... second input
    //! @return ... result
    T Evaluate(const T& rInput1, const T& rInput2)const
    {
        if (rInput1<rInput2)
            return rInput1;
        return rInput2;
    }
};

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all Operators in NuTo
template <class T>
class DOperatorMax : public DyadicOperator<T>
{
public:
    DOperatorMax<T>() :  DyadicOperator<T>()
    {}

    //! @brief ... operator with two inputs and one output
    //! @param rInput1 ... first input
    //! @param rInput2 ... second input
    //! @return ... result
    T Evaluate(const T& rInput1, const T& rInput2)const
    {
        if (rInput1>rInput2)
            return rInput1;
        return rInput2;
    }
};

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all Operators in NuTo
template <class T>
class DOperatorAdd : public DyadicOperator<T>
{
public:
    DOperatorAdd<T>() :  DyadicOperator<T>()
    {}

    //! @brief ... operator with two inputs and one output
    //! @param rInput1 ... first input
    //! @param rInput2 ... second input
    //! @return ... result
    T Evaluate(const T& rInput1, const T& rInput2)const
    {
        return rInput1 + rInput2;
    }
};

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all Operators in NuTo
template <class T>
class DOperatorSub : public DyadicOperator<T>
{
public:
    DOperatorSub<T>() :  DyadicOperator<T>()
    {}

    //! @brief ... operator with two inputs and one output
    //! @param rInput1 ... first input
    //! @param rInput2 ... second input
    //! @return ... result
    T Evaluate(const T& rInput1, const T& rInput2)const
    {
        return rInput1 - rInput2;
    }
};

//! @author Jörg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all Operators in NuTo
template <class T>
class DOperatorMul : public DyadicOperator<T>
{
public:
    DOperatorMul<T>() :  DyadicOperator<T>()
    {}

    //! @brief ... operator with two inputs and one output
    //! @param rInput1 ... first input
    //! @param rInput2 ... second input
    //! @return ... result
    T Evaluate(const T& rInput1, const T& rInput2)const
    {
        return rInput1 * rInput2;
    }
};

} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::MonadicOperator<int>)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::MonadicOperator<double>)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::DyadicOperator<int>)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::DyadicOperator<double>)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
#endif // OPERATOR_H
