#ifndef METAMODEL_H
#define METAMODEL_H

#include <string>
#include <time.h>

#include "nuto/metamodel/MetamodelException.h"
#include "nuto/base/NuToObject.h"
#include "nuto/metamodel/SupportPoints.h"

extern "C" { 
#include <dSFMT.h>
}

namespace NuTo
{
template<class T> class FullMatrix;

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all metamodels in NuTo
class Metamodel : public virtual NuToObject
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    Metamodel() : NuToObject()
    {
        // init random number generator with milliseconds from ..
        dsfmt_init_gen_rand(&mRandomNumberGenerator, time (NULL));
    }
    
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {    
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
           & BOOST_SERIALIZATION_NVP(mSupportPoints);
    }
#endif  // ENABLE_SERIALIZATION

    void Build();
    virtual void BuildDerived()=0;
    
    void AppendMinMaxTransformationInput(double min, double max);
    void AppendMinMaxTransformationInput(int coordinate, double min, double max);
    void AppendMinMaxTransformationOutput(double min, double max);
    void AppendMinMaxTransformationOutput(int coordinate, double min, double max);

    //! @brief ... append for all coordinates of the input points a transformation to zero mean and unit variance to the list of transformations.
    void AppendZeroMeanUnitVarianceTransformationInput();
    
    //! @brief ... append for a specific coordinate of the input points a transformation to zero mean and unit variance to the list of transformations.
    //! @param rCoordinate ... coordinate for which the transformation is performed
    void AppendZeroMeanUnitVarianceTransformationInput(int rCoordinate);

    //! @brief ... append for all coordinates of the output points a transformation to zero mean and unit variance to the list of transformations.
    void AppendZeroMeanUnitVarianceTransformationOutput();
    
    //! @brief ... append for a specific coordinate of the output points a transformation to zero mean and unit variance to the list of transformations.
    //! @param rCoordinate ... coordinate for which the transformation is performed
    void AppendZeroMeanUnitVarianceTransformationOutput(int rCoordinate);

    FullMatrix<double> GetOriginalSupportPointsInput()const;
    FullMatrix<double> GetOriginalSupportPointsOutput()const;
    FullMatrix<double> GetTransformedSupportPointsInput()const;
    FullMatrix<double> GetTransformedSupportPointsOutput()const;
    void SetSupportPoints(int rDimInput, int rDimOutput, FullMatrix<double> rInputCoordinates, FullMatrix<double> rOutputCoordinates);
    void BuildTransformation();
    void InitRandomNumberGenerator(int rSeed);
    double RandomDouble();
    virtual void Info()const;
    void Solve(const FullMatrix<double>& rInputCoordinates, FullMatrix<double>& rOutputCoordinates)const;
    void SolveConfidenceInterval(const FullMatrix<double>& rInputCoordinates, NuTo::FullMatrix<double>& rOutputCoordinates, 
                                                    NuTo::FullMatrix<double>& rOutputCoordinatesMin, NuTo::FullMatrix<double>& rOutputCoordinatesMax)const;
                                                    virtual void SolveTransformed(const FullMatrix<double>& rInputCoordinates, NuTo::FullMatrix<double>& rOutputCoordinates)const=0;    
    virtual void SolveConfidenceIntervalTransformed(const FullMatrix<double>& rInputCoordinates, NuTo::FullMatrix<double>& rOutputCoordinates, 
                                            NuTo::FullMatrix<double>& rOutputCoordinatesMin, NuTo::FullMatrix<double>& rOutputCoordinatesMax)const
    {
        throw MetamodelException("Metamodel::SolveConfidenceIntervalTransformed - not implemented for this kind of metamodel.");
    }
protected:

    NuTo::SupportPoints mSupportPoints;
    dsfmt_t mRandomNumberGenerator;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
#include <boost/serialization/assume_abstract.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Metamodel)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
#endif // METAMODEL_H
