// $Id$

#ifndef NystroemBase_H
#define NystroemBase_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for explicit time integration using Nyström methods
//! see for example Hairer, Lubich and Wanner ("Geometrical numerical integration")
class NystroemBase : public TimeIntegrationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NystroemBase(StructureBase& rStructure);

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif// SWIG
#endif // ENABLE_SERIALIZATION

    //! @brief perform the time integration
    //! @param rStructure ... structure
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(StructureBase& rStructure, double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief ... Return time step
    double GetTimeStep()const
    {
    	return mTimeStep;
    }

    //! @brief ... Return time step
    void SetTimeStep(double rTimeStep)
    {
    	mTimeStep = rTimeStep;
    }

    //! @brief ... return number of intermediate stages
    virtual int GetNumStages()const=0;

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    virtual double GetStageTimeFactor(int rStage)const=0;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    virtual void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const=0;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    virtual double GetStageWeights1(int rStage)const=0;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    //! for standard Runge Kutta-Methods, this is identical
    //! for Nyström methods, this migh be different
    virtual double GetStageWeights2(int rStage)const=0;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    virtual bool HasTimeChanged(int rStage)const=0;

protected:
	//time step for the time integration, be careful not to make it smaller than the critical time step
    double mTimeStep;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::NystroemBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



#endif // NystroemBase_H
