#pragma once

#include "metamodel/Metamodel.h"
#include "optimize/CallbackHandler.h"

#include "base/Exception.h"
#include <Eigen/Core>
#include <string>

namespace NuTo
{

class TransferFunction;

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief Standard abstract class for all metamodels in NuTo
class NeuralNetwork : public Metamodel, public CallbackHandler
{
public:
    enum eTransferFunctions
    {
        Empty = 0,
        HardLim,
        HardLims,
        PureLin,
        SatLin,
        SatLins,
        LogSig,
        TanSig,
        PosLin
    };

    //! @brief constructor
    //! @param rvNumNeurons ... number of neurons in each hidden layer
    NeuralNetwork(std::vector<int> rvNumNeurons);

    void BuildDerived() override;
    void SetTransferFunction(int rLayer, eTransferFunctions rTransferFunction);

    double Objective() const override;
    void Gradient(Eigen::Ref<Eigen::MatrixXd> rGradient) const override;
    void Hessian(Eigen::MatrixXd& rDiagHessian) const override;
    void HessianDiag(Eigen::MatrixXd& rDiagHessian) const;
    void HessianFull(Eigen::Ref<Eigen::MatrixXd> rDiagHessian) const;

    Eigen::MatrixXd GetParameters() const;
    void SetParameters(const Eigen::MatrixXd& Parameters) override;
    virtual void Info() const override;

    //! @brief ... get the inverse noise covariance matrix
    //! @param rInverseCovariance ... inverse noise covariance matrix
    void GetInverseNoiseCovarianceMatrixTransformed(Eigen::MatrixXd& rInverseCovariance) const;

    //! @brief ... get the noise covariance matrix
    //! @param rCovariance ... noise covariance matrix
    void GetNoiseCovarianceMatrixTransformed(Eigen::MatrixXd& rCovariance) const;

    //! @brief ... get noise correlation matrix
    //! @param rNoiseCorrelation ... noise correlation matrix
    void GetNoiseCorrelationMatrix(Eigen::MatrixXd& rNoiseCorrelation) const;

    //! @brief ... get precision (alpha) parameters
    //! @param rPrecisionParameters ... precision parameters
    void GetPrecisionParametersTransformed(Eigen::MatrixXd& rPrecisionParameters) const;

    int GetNumParameters() const
    {
        return (mNumWeights + mNumBiases);
    }

    int GetNumNeurons() const
    {
        int numNeurons(0);
        for (int currentLayer = 0; currentLayer < mNumLayers + 1; currentLayer++)
            numNeurons += mvNumNeurons[currentLayer];
        return numNeurons;
    }

    void SetInitAlpha(const double rInitAlpha)
    {
        if (rInitAlpha < 0)
            throw Exception("NuTo::NeuralNetwork::SetInitAlpha - Init alpha must be non negative.");
        mInitAlpha = rInitAlpha;
    }

    void SetAccuracyGradient(const double rAccuracyGradient)
    {
        mAccuracyGradient = rAccuracyGradient;
    }

    void SetMinObjective(const double rMinObjective)
    {
        mMinObjective = rMinObjective;
    }

    void SetMinDeltaObjectiveBetweenRestarts(const double rMinDeltaObjectiveBetweenRestarts)
    {
        mMinDeltaObjectiveBetweenRestarts = rMinDeltaObjectiveBetweenRestarts;
    }
    void SetMinDeltaObjectiveBayesianIteration(const double rMinDeltaObjectiveBayesianIteration)
    {
        mMinDeltaObjectiveBayesianIteration = rMinDeltaObjectiveBayesianIteration;
    }

    void SetMaxFunctionCalls(const int rMaxFunctionCalls)
    {
        mMaxFunctionCalls = rMaxFunctionCalls;
    }

    void SetShowSteps(const int rShowSteps)
    {
        mShowSteps = rShowSteps;
    }

    void SetMaxBayesianIterations(const int rMaxBayesianIterations)
    {
        mMaxBayesianIterations = rMaxBayesianIterations;
    }


    inline void UseDiagHessian()
    {
        mUseDiagHessian = true;
    }

    inline void UseFullHessian()
    {
        mUseDiagHessian = false;
    }

    inline void SetBayesianTraining()
    {
        mBayesian = true;
    }

    inline void UnsetBayesianTraining()
    {
        mBayesian = false;
    }

    void Jacobian(Eigen::MatrixXd& rJacobian, std::vector<double>& pA, std::vector<double>& pO,
                  Eigen::MatrixXd& pM) const;

    void SolveTransformed(const Eigen::MatrixXd& rInputCoordinates, Eigen::MatrixXd& rOutputCoordinates) const override;
    void SolveConfidenceIntervalTransformed(const Eigen::MatrixXd& rInputCoordinates,
                                            Eigen::MatrixXd& rOutputCoordinates, Eigen::MatrixXd& rOutputCoordinatesMin,
                                            Eigen::MatrixXd& rOutputCoordinatesMax) const override;

protected:
    void ForwardPropagateInput(std::vector<double>& pA, std::vector<double>& pO) const;
    void GetAlphas(Eigen::VectorXd& rAlpha) const; // calculate for each free parameter the corresponding alpha
    void GetPosInAlphaVector(std::vector<int>& rPosInAlphaVector) const;
    void GetRefsPerAlpha(std::vector<int>& rRefsPerAlpha) const;
    void dLnDetA_dBeta(Eigen::MatrixXd& rHessianInv, Eigen::MatrixXd& rResult) const;

    //! @brief number of layers (excluding input and including output)
    int mNumLayers;
    //! @brief decide, if a Bayesian training procedure should be applied
    bool mBayesian;
    //! @brief decide, that when calling the hessian matrix only the diagonal values are calculated
    bool mUseDiagHessian;
    //! @brief transfer functions for each hidden layer including output
    std::vector<TransferFunction*> mvTransferFunction;
    //! @brief number of neurons in each hidden layer + input and output layer for compatibility
    std::vector<int> mvNumNeurons;
    //! @brief vector of weights
    /*!
        size is mNumLayer, for each connection between layers the corresponding weights
        first all the connections of a single neuron in the current layer to all the neurons in previous layer,
        and then next neuron in current layer
    */
    std::vector<std::vector<double>> mvWeights;
    //! @brief total number of weights
    int mNumWeights;
    //! @brief total number of biases
    int mNumBiases;
    //! @brief vector of biases
    /*!
        size is mNumLayer, for each layer the corresponding biases
    */
    std::vector<std::vector<double>> mvBias;
    //! @brief inverse convariance between the error terms, size is DimOutput*DimOutput
    Eigen::MatrixXd mvCovarianceInv;
    //! @brief accuracy of weights
    /*!
        size is DimInput+DimOutput+DimHiddenLayer
        weights inputs, bias first hiddenlayer, weights/biases second,..,weights and biases for each output
     */
    std::vector<double> mvAlpha;
    //! @brief Initialization constant for the alphas
    double mInitAlpha;
    //! @brief Stopping criteria for the optimizer (norm of the gradient)
    double mAccuracyGradient;
    //! @brief Stopping criteria for the optimizer (smallest objective)
    double mMinObjective;
    //! @brief Stopping criteria for the optimizer (minimum delta objective between restarts of the CG algorithm)
    double mMinDeltaObjectiveBetweenRestarts;
    //! @brief Stopping criteria for the Bayesian iteration (relative delta objective between two consecutive
    //! iterations)
    double mMinDeltaObjectiveBayesianIteration;
    //! @brief Stopping criteria for the optimizer (maximum number of function calls)
    int mMaxFunctionCalls;
    //! @brief number of steps to show during the optimization phase
    int mShowSteps;
    //! @brief maximum number of iterations for the optimization of the hyperparameters
    int mMaxBayesianIterations;

    //! @brief constructor required for serialization
    NeuralNetwork();
};
} // namespace NuTo
