#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "metamodel/Metamodel.h"
#include "optimize/CallbackHandler.h"

#include "metamodel/MetamodelException.h"
#include <eigen3/Eigen/Core>
#include <string>

namespace NuTo
{

class TransferFunction;

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief Standard abstract class for all metamodels in NuTo
class NeuralNetwork : public virtual Metamodel, public virtual CallbackHandler
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
enum eTransferFunctions {
       Empty=0,
	   HardLim,
	   HardLims,
	   PureLin,
	   SatLin,
	   SatLins,
	   LogSig,
	   TanSig,
	   PosLin};

	//! @brief constructor
	//! @param rvNumNeurons ... number of neurons in each hidden layer
    NeuralNetwork (std::vector<int> rvNumNeurons);

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class, this is the load routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    //! @brief serializes the class, this is the save routine
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
#endif  // ENABLE_SERIALIZATION

	void BuildDerived();
	void SetTransferFunction(int rLayer, eTransferFunctions rTransferFunction);
	
	double Objective()const;
	void Gradient(Eigen::MatrixXd& rGradient)const;
    void Hessian(Eigen::MatrixXd&  rDiagHessian)const;
    void HessianDiag(Eigen::MatrixXd&  rDiagHessian)const;
    void HessianFull(Eigen::MatrixXd&  rDiagHessian)const;
    
    Eigen::MatrixXd GetParameters() const;
    void SetParameters(const Eigen::MatrixXd& Parameters);
    virtual void Info()const;
    
    //! @brief ... get the inverse noise covariance matrix
    //! @param rInverseCovariance ... inverse noise covariance matrix
    void GetInverseNoiseCovarianceMatrixTransformed(Eigen::MatrixXd& rInverseCovariance)const;
    
    //! @brief ... get the noise covariance matrix
    //! @param rCovariance ... noise covariance matrix
    void GetNoiseCovarianceMatrixTransformed(Eigen::MatrixXd& rCovariance)const;
    
    //! @brief ... get noise correlation matrix
    //! @param rNoiseCorrelation ... noise correlation matrix
    void GetNoiseCorrelationMatrix(Eigen::MatrixXd& rNoiseCorrelation)const;
    
    //! @brief ... get precision (alpha) parameters
    //! @param rPrecisionParameters ... precision parameters
    void GetPrecisionParametersTransformed(Eigen::MatrixXd& rPrecisionParameters)const;
    
    int GetNumParameters()const
	{
	    return (mNumWeights + mNumBiases);
	}

	int GetNumNeurons()const
	{
    	int numNeurons(0);
		for (int currentLayer=0; currentLayer<mNumLayers+1; currentLayer++)
        	numNeurons+=mvNumNeurons[currentLayer];
		return numNeurons;
	}
    
    void SetInitAlpha(const double rInitAlpha)
    {
        if (rInitAlpha<0)
            throw MetamodelException("NuTo::NeuralNetwork::SetInitAlpha - Init alpha must be non negative."); 
        mInitAlpha=rInitAlpha; 
    }

    void SetAccuracyGradient(const double rAccuracyGradient)
    {
        mAccuracyGradient=rAccuracyGradient; 
        
    }

    void SetMinObjective(const double rMinObjective)
    {
        mMinObjective=rMinObjective; 
        
    }

    void SetMinDeltaObjectiveBetweenRestarts(const double rMinDeltaObjectiveBetweenRestarts)
    {
        mMinDeltaObjectiveBetweenRestarts=rMinDeltaObjectiveBetweenRestarts; 
        
    }
    void SetMinDeltaObjectiveBayesianIteration(const double rMinDeltaObjectiveBayesianIteration)
    {
        mMinDeltaObjectiveBayesianIteration=rMinDeltaObjectiveBayesianIteration; 
        
    }
    
    void SetMaxFunctionCalls(const int rMaxFunctionCalls)
    {
        mMaxFunctionCalls=rMaxFunctionCalls; 
        
    }

    void SetShowSteps(const int rShowSteps)
    {
        mShowSteps=rShowSteps; 
        
    }

    void SetMaxBayesianIterations(const int rMaxBayesianIterations)
    {
        mMaxBayesianIterations=rMaxBayesianIterations; 
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
    
#ifndef SWIG
    void Jacobian(Eigen::MatrixXd& rJacobian, std::vector<double>& pA, std::vector<double>& pO, Eigen::MatrixXd& pM)const;
#endif
    
    void SolveTransformed(const Eigen::MatrixXd& rInputCoordinates, Eigen::MatrixXd& rOutputCoordinates) const;
    void SolveConfidenceIntervalTransformed(const Eigen::MatrixXd& rInputCoordinates, Eigen::MatrixXd& rOutputCoordinates,
            Eigen::MatrixXd& rOutputCoordinatesMin, Eigen::MatrixXd& rOutputCoordinatesMax) const;
#ifdef ENABLE_SERIALIZATION
    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    //! @brief ... save the object to a file
    virtual void Restore (const std::string &filename, std::string rType );

	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
	virtual void Save (const std::string &filename, std::string rType )const;
#endif // ENABLE_SERIALIZATION

	//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const;

protected:                                                                 
    void ForwardPropagateInput(std::vector<double> &pA, std::vector<double> &pO)const;
    void GetAlphas (Eigen::VectorXd& rAlpha)const;               //calculate for each free parameter the corresponding alpha
    void GetPosInAlphaVector (std::vector<int>& rPosInAlphaVector)const;
    void GetRefsPerAlpha (std::vector<int>& rRefsPerAlpha)const;
    void dLnDetA_dBeta(Eigen::MatrixXd& rHessianInv, Eigen::MatrixXd& rResult)const;
    
    //! @brief number of layers (excluding input and including output)
    int mNumLayers;
    //! @brief decide, if a Bayesian training procedure should be applied
    bool mBayesian;
    //! @brief decide, that when calling the hessian matrix only the diagonal values are calculated
    bool mUseDiagHessian;
    //! @brief transfer functions for each hidden layer including output
    std::vector<TransferFunction*>    mvTransferFunction;
    //! @brief number of neurons in each hidden layer + input and output layer for compatibility
    std::vector<int>                  mvNumNeurons;
    //! @brief vector of weights
    /*! 
        size is mNumLayer, for each connection between layers the corresponding weights 
        first all the connections of a single neuron in the current layer to all the neurons in previous layer, 
        and then next neuron in current layer 
    */
    std::vector<std::vector<double> > mvWeights;
    //! @brief total number of weights
	int                               mNumWeights;
    //! @brief total number of biases
	int                               mNumBiases;
    //! @brief vector of biases
    /*!
        size is mNumLayer, for each layer the corresponding biases
    */
    std::vector<std::vector<double> > mvBias;
    //! @brief inverse convariance between the error terms, size is DimOutput*DimOutput
    Eigen::MatrixXd                   mvCovarianceInv;
    //! @brief accuracy of weights
    /*!
        size is DimInput+DimOutput+DimHiddenLayer
        weights inputs, bias first hiddenlayer, weights/biases second,..,weights and biases for each output
     */
    std::vector<double>               mvAlpha;
    //! @brief Initialization constant for the alphas
    double                            mInitAlpha;
    //! @brief Stopping criteria for the optimizer (norm of the gradient)
    double                            mAccuracyGradient;
    //! @brief Stopping criteria for the optimizer (smallest objective)
    double                            mMinObjective;
    //! @brief Stopping criteria for the optimizer (minimum delta objective between restarts of the CG algorithm)
    double                            mMinDeltaObjectiveBetweenRestarts;
    //! @brief Stopping criteria for the Bayesian iteration (relative delta objective between two consecutive iterations)
    double                            mMinDeltaObjectiveBayesianIteration;
    //! @brief Stopping criteria for the optimizer (maximum number of function calls)
    int                               mMaxFunctionCalls;
    //! @brief number of steps to show during the optimization phase
    int                               mShowSteps;
    //! @brief maximum number of iterations for the optimization of the hyperparameters
    int                               mMaxBayesianIterations;

    //! @brief constructor required for serialization
    NeuralNetwork ();
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::NeuralNetwork)
#endif //SWIG
#endif  // ENABLE_SERIALIZATION
