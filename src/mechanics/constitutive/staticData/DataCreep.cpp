#include "mechanics/constitutive/staticData/DataCreep.h"

#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"


using namespace NuTo::Constitutive::StaticData;


template <typename TStream>
void NuTo::Constitutive::StaticData::DataCreep::SerializeDataCreep(TStream& rStream)
{
}


template void DataCreep::SerializeDataCreep<NuTo::SerializeStreamIn>(SerializeStreamIn& rStream);
template void DataCreep::SerializeDataCreep<NuTo::SerializeStreamOut>(SerializeStreamOut& rStream);

void DataCreep::ProceedToNextTimestep()
{
    mPreviousTime = mCurrentTime;
}
