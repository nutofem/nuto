
/*******************************************************************************
 Bauhaus-Universitaet Weimar
 Author: Joerg F. Unger,  September 2009
*******************************************************************************/

#pragma once

namespace NuTo
{

class TransferFunction
{

public:
    typedef enum {
        EMPTY_TRANSFER_FUNCTION,
        HardLim,
        HardLimS,
        PURELIN,
        SATLIN,
        SATLINS,
        LOGSIG,
        TANSIG,
        POSLIN
    } eTransferFunction;

    TransferFunction()
    {
    }

    virtual ~TransferFunction()
    {
    }

    virtual double evaluate(double x) = 0;
    virtual TransferFunction* clone() const = 0;
    virtual double derivative(double x) = 0;
    virtual double second_derivative(double x) = 0;
    virtual void info() const = 0;

    virtual eTransferFunction get_enum() const = 0;
};

/*******************************************************************************
 none
 *******************************************************************************/
class EmptyTransferFunction : public TransferFunction
{
public:
    EmptyTransferFunction()
    {
    }

    ~EmptyTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return EMPTY_TRANSFER_FUNCTION;
    }
};

/*******************************************************************************
 hard limit function
 *******************************************************************************/
class HardLimTransferFunction : public TransferFunction
{

public:
    HardLimTransferFunction()
    {
    }

    ~HardLimTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return HardLim;
    }
};

/*******************************************************************************
 symmetric hard limit function
 *******************************************************************************/
class HardLimsTransferFunction : public TransferFunction
{

public:
    HardLimsTransferFunction()
    {
    }

    ~HardLimsTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return HardLimS;
    }
};

/*******************************************************************************
 linear limit function
 *******************************************************************************/
class PureLinTransferFunction : public TransferFunction
{

public:
    PureLinTransferFunction()
    {
    }

    ~PureLinTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return PURELIN;
    }
};

/*******************************************************************************
 saturating linear limit function
 *******************************************************************************/
class SatLinTransferFunction : public TransferFunction
{

public:
    SatLinTransferFunction()
    {
    }

    ~SatLinTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return SATLIN;
    }
};

/*******************************************************************************
 symmetric saturating linear limit function
 *******************************************************************************/
class SatLinsTransferFunction : public TransferFunction
{

public:
    SatLinsTransferFunction()
    {
    }

    ~SatLinsTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return SATLINS;
    }
};

/*******************************************************************************
 log sigmoid function
 *******************************************************************************/
class LogSigTransferFunction : public TransferFunction
{
public:
    LogSigTransferFunction()
    {
    }

    ~LogSigTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return LOGSIG;
    }
};

/*******************************************************************************
 hyperbolic tangent sigmoid
 *******************************************************************************/
class TanSigTransferFunction : public TransferFunction
{
public:
    TanSigTransferFunction()
    {
    }

    ~TanSigTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return TANSIG;
    }
};

/*******************************************************************************
 positive linear
 *******************************************************************************/
class PosLinTransferFunction : public TransferFunction
{
public:
    PosLinTransferFunction()
    {
    }

    ~PosLinTransferFunction()
    {
    }

    double evaluate(double x) override;
    TransferFunction* clone() const override;
    double derivative(double x) override;
    double second_derivative(double x) override;
    void info() const override;

    TransferFunction::eTransferFunction get_enum() const override
    {
        return POSLIN;
    }
};

} // namespace surrogate
