#include "nuto/mechanics/constitutive/staticData/Component.h"
#include "nuto/base/serializeStream/SerializeStreamIn_Def.h"
#include "nuto/base/serializeStream/SerializeStreamOut_Def.h"

namespace NuTo
{
namespace Constitutive
{
namespace StaticData
{


SerializeStreamOut& operator<<(SerializeStreamOut& rStream, Component& rData)
{
    rData.WriteComponent(rStream);
    return rStream;
}

SerializeStreamIn& operator>>(SerializeStreamIn& rStream,Component& rData)
{
    rData.ReadComponent(rStream);
    return rStream;
}
}
}
}
