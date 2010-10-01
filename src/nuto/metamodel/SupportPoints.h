// $Id$

/*******************************************************************************
Bauhaus-Universität Weimar
Author: Joerg F. Unger,  Septermber 2009
*******************************************************************************/


#ifndef SUPPORTPOINTS_H
#define SUPPORTPOINTS_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION
#include <boost/ptr_container/ptr_list.hpp>

#include "nuto/metamodel/Transformation.h"
#include "nuto/math/FullMatrix.h"

namespace NuTo
{
template<class T> class FullMatrix;
//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief stores the support points
class SupportPoints
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
	//! @brief constructor
    SupportPoints();
    //! @brief destructor
    ~SupportPoints();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief clear support points
    void Clear();

    //! @brief info about support points
    void Info()const;

    //! @brief get number of support points
    inline int GetNumSupportPoints()const
	{
        return mSPOrigInput.GetNumColumns();
	}

    //! @brief get input dimension of support points (input)
    inline int GetDimInput()const                              
	{
        return mSPOrigInput.GetNumRows();
	}

    //! @brief get output dimension of support points (input)
    inline int GetDimOutput()const
	{
        return mSPOrigOutput.GetNumRows();
	}

    //! @brief returns the input of the support points in a matrix
    inline const FullMatrix<double>& GetOrigSupportPointsInput()const
	{
	    return mSPOrigInput;
	}
    
    //! @brief returns the output of the support points in a matrix
    inline const FullMatrix<double>& GetOrigSupportPointsOutput()const
	{
	    return mSPOrigOutput;
	}
    
    //! @brief returns the input of the transformed support points in a matrix
    inline const FullMatrix<double>& GetTransformedSupportPointsInput()const
	{
	    return mSPTransInput;
	}
    
    //! @brief returns the input of the transformed support points in a matrix
    inline const FullMatrix<double>& GetTransformedSupportPointsOutput()const
	{
	    return mSPTransOutput;
	}
    
    //! @brief append a transformation for input points
    //! @brief in general, the order is orig_input -> outputtrans1 -> outputtrans2 -> trans_input
    void AppendTransformationInput(Transformation* rTransformation);

    //! @brief append a transformation for output points
    //! @brief in general, order is orig_output -> outputtrans1 -> outputtrans2 -> trans_output
    void AppendTransformationOutput(Transformation* rTransformation);

    //! @brief checks, if the transformation has been build
	bool IsTransformationBuild()const
	{
	    return mTransformationBuild;
	}

    //! @brief build transformation for support points
    void BuildTransformation();

    //! @brief set support points
    void SetSupportPoints(const FullMatrix<double>& rSPOrigInput, const FullMatrix<double>& rSPOrigOutput);

    //! @brief return the Weight vector as Matrix (use only from Python level, since everything is copied)
    inline const FullMatrix<double> GetWeights(int rSample)const
    {
        return FullMatrix<double>(mWeight.size(),1,mWeight);
    }

    //! @brief Removes all applied Transformations for the input and the output points
    void ClearTransformations();
    
#ifndef SWIG
    inline double GetWeight(int rSample)const
	{
        return mWeight[rSample];
	}

    //! @brief perform forward transformation for inputs (from orig to transformed)
    //! @brief in general, the forward direction is orig_input -> inputtrans1 ->inputtrans2 -> metamodel -> outputtrans2 -> outputtrans1 -> orig_output
    void TransformForwardInput(FullMatrix<double>& rCoordinates)const;
    
    //! @brief perform backward transformation for inputs (from transformed to orig
    void TransformBackwardInput(FullMatrix<double>& rCoordinates)const;
    
    //! @brief perform forward transformation for outputs (from transformed to orig)
    //! @brief attention, this is exactly the backwards order, since transformations for outputs are given in revers order
    void TransformForwardOutput(FullMatrix<double>& rCoordinates)const;
    
    //! @brief perform backward transformation for outputs
    void TransformBackwardOutput(FullMatrix<double>& rCoordinates)const;
    
#endif

private:
    FullMatrix<double>  mSPOrigInput;      //!< original inputs, each sample after another
    FullMatrix<double>  mSPOrigOutput;     //!< original outputs, each sample after another
    FullMatrix<double>  mSPTransInput;     //!< transformed inputs, each sample after another
    FullMatrix<double>  mSPTransOutput;    //!< transformed outputs, each sample after another

    std::vector<double> mWeight;    //!< weight of each support point
    
    boost::ptr_list<Transformation> mlTransformationInput;   //!< list of Transformations for inputs
    boost::ptr_list<Transformation> mlTransformationOutput;  //!< list of Transformations for outputs
	
	bool mTransformationBuild;
    


};
} // namespace nuto
#endif /* SUPPORTPOINTS_H */
