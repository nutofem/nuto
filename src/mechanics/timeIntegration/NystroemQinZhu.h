// $Id$

#pragma once

#include "mechanics/timeIntegration/NystroemBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, NU
//! @date November 2013
//! @brief ... class for explicit time integration using a fourth order Nystroem method according to
//! "Modeling the scalar wave equation with Nyström methods" by Jing-Bo Chen, Geophysics, 2006
//! "High-order time discretizations in seismic modeling", Jing-Bo Chen 2007
class NystroemQinZhu : public NystroemBase
{

public:

    //! @brief constructor
    NystroemQinZhu(StructureBase* rStructure);

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep()const override
    {
    	return true;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep()const override;

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const override;

    //! @brief ... return number of intermediate stages
    int GetNumStages()const override
    {
    	return 3;
    }

    //! @brief ... return delta time factor of intermediate stages (c in Butcher tableau)
    double GetStageTimeFactor(int rStage)const override;

    //! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
    void GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const override;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    double GetStageWeights1(int rStage)const override;

    //! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
    double GetStageWeights2(int rStage)const override;

    //! @brief ... return, if time (e.g. for the calculation of external loads) has changed
    bool HasTimeChanged(int rStage)const override;

protected:
    //empty private construct required for serialization
};
} //namespace NuTo



