// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/timeIntegration/TimeIntegrationBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (Newmark, but you can use it for statics as well with setting the flag isDynamic to false)
class RungeKuttaBase : public TimeIntegrationBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    RungeKuttaBase(StructureBase* rStructure);

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
    //! @param rTimeDelta ... length of the simulation
    NuTo::eError Solve(double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

    //! @brief ... return number of intermediate stages
    virtual int GetNumStages()const=0;

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    virtual double GetStageTimeFactor(int rStage)const=0;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    virtual void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const=0;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    virtual double GetStageWeights(int rStage)const=0;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    virtual bool HasTimeChanged(int rStage)const=0;

protected:
    //empty private construct required for serialization
#ifdef ENABLE_SERIALIZATION
    RungeKuttaBase(){};
#endif  // ENABLE_SERIALIZATION
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::RungeKuttaBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION



