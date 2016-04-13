/*
 * ImplEx.cpp
 *
 *  Created on: 17 Mar 2016
 *      Author: ttitsche
 */

#include "nuto/mechanics/timeIntegration/ImplEx.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/structures/StructureOutputDummy.h"


#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

NuTo::ImplEx::ImplEx(StructureBase* rStructure) : ImplicitExplicitBase(rStructure)
{
    mExtrapolationErrorThreshold = -1.;
}

void NuTo::ImplEx::ExtrapolateStaticData(const ConstitutiveTimeStep& rTimeStep)
{
    // Setup input list with the timeStep and the option to use EULER_FORWARD
    ConstitutiveInputMap input;
    ConstitutiveCalculateStaticData calculateStaticDataEulerForward(CalculateStaticData::EULER_FORWARD);
    ConstitutiveTimeStep nonconstTimeStepCopy(rTimeStep);

    input[Constitutive::Input::CALCULATE_STATIC_DATA] = &calculateStaticDataEulerForward;
    input[Constitutive::Input::TIME_STEP] = &nonconstTimeStepCopy;


    // Setup output map to store the static data
    std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> evalUpdateStaticData;
    StructureOutputDummy dummy;
    evalUpdateStaticData[StructureEnum::UPDATE_STATIC_DATA] = &dummy;

    mStructure->Evaluate(input, evalUpdateStaticData);
}

bool NuTo::ImplEx::CheckExtrapolationAndAdjustTimeStep()
{
    if (mExtrapolationErrorThreshold == -1.)
        throw MechanicsException(__PRETTY_FUNCTION__, "Define the extrapolation threshold first.");

    double extrapolationError = mStructure->ElementTotalGetStaticDataExtrapolationError();
    double increaseTimeStepThreshold = mExtrapolationErrorThreshold / 10;

    bool extrapolationIsOK = extrapolationError < mExtrapolationErrorThreshold;

    std::cout << "relative Extrapolation error (>1 is bad): " << extrapolationError / mExtrapolationErrorThreshold << std::endl;

    if (not extrapolationIsOK && not mForceAcceptOfNextSolution)
    {
        mTimeStep *= .1;
        mStructure->GetLogger() << "[" << __FUNCTION__ << "] ### decreasing time step to " << mTimeStep << "\n";
        mForceAcceptOfNextSolution = true;
//        std::cin.ignore();
        return false;
    }

    // accept solution, either the extrapolation is OK or this solution is forced to be accepted (somehow continue the time stepping)

    mForceAcceptOfNextSolution = false;

    if (extrapolationError < increaseTimeStepThreshold)
    {
        mTimeStep *= 1.3;
        mStructure->GetLogger() << "[" << __FUNCTION__ << "] increasing time step to " << mTimeStep << "\n";
    }

    return true;
}



#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ImplEx::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ImplEx::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ImplEx::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of ImplEx" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase);

    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of ImplEx" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION
