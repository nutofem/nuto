// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/RungeKuttaBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date February 2012
//! @brief ... standard class for implicit timeintegration (Newmark, but you can use it for statics as well with setting
//! the flag isDynamic to false)
class RungeKuttaDormandPrince : public RungeKuttaBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    RungeKuttaDormandPrince(StructureBase* rStructure);

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep() const
    {
        return true;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep() const;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in
    //! the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId() const;

    // returns the accelerations
    void f(NuTo::StructureBase* mStructure, NuTo::SparseDirectSolverMUMPS& mySolver,
           const NuTo::StructureOutputBlockVector& extLoad, const NuTo::StructureOutputBlockVector& dof_dt0,
           const NuTo::StructureOutputBlockVector& dof_dt1, double factor,
           const NuTo::StructureOutputBlockVector& rAcceleration0, const NuTo::StructureOutputBlockVector& rVelocity0,
           NuTo::StructureOutputBlockVector& rAcceleration, NuTo::StructureOutputBlockVector& rVelocity, int stage);

    NuTo::eError DP_DoStep(double rTimeDelta, int rLoadCase, double curTime, NuTo::SparseDirectSolverMUMPS& mySolver,
                           std::vector<StructureOutputBlockVector>& kAcc, std::vector<StructureOutputBlockVector>& kVel,
                           StructureOutputBlockVector& extLoad, NuTo::StructureOutputBlockVector& dof_dt0,
                           NuTo::StructureOutputBlockVector& dof_dt1);

    //! @brief ... return number of intermediate stages
    int GetNumStages() const
    {
        return 7;
    }

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    double GetStageTimeFactor(int rStage) const;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage) const;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    double GetStageWeights(int rStage) const;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    bool HasTimeChanged(int rStage) const;

protected:
// empty private construct required for serialization
#ifdef ENABLE_SERIALIZATION
    RungeKuttaDormandPrince(){};
#endif // ENABLE_SERIALIZATION
};
} // namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::RungeKuttaDormandPrince)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
