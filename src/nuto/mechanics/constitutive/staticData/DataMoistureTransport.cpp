#include "nuto/mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"
#include "nuto/base/serializeStream/SerializeStreamOut.h"

using namespace NuTo::Constitutive::StaticData;

NuTo::FullVector<double,Eigen::Dynamic> DataMoistureTransport::GetCurrentSorptionCoeff () const
{
    return mCurrentSorptionCoeff;
}


double DataMoistureTransport::GetLastRelHumValue() const
{
    return mLastRelHumValue;
}


NuTo::FullVector<double,Eigen::Dynamic> DataMoistureTransport::GetLastSorptionCoeff() const
{
    return mLastSorptionCoeff;
}


bool DataMoistureTransport::IsDesorption() const
{
    return mSorptionHistoryDesorption;
}


void DataMoistureTransport::SetCurrentSorptionCoeff(FullVector<double,Eigen::Dynamic> rCurrentSorptionCoeff)
{
    switch (rCurrentSorptionCoeff.GetNumRows())
    {
    case 3:
    {
        mCurrentSorptionCoeff = rCurrentSorptionCoeff;
        break;
    }
    case 4:
    {
        mCurrentSorptionCoeff.Resize(3);
        for(int i=0; i<3; i++)
        {
            mCurrentSorptionCoeff(i) = rCurrentSorptionCoeff(i+1);
        }
        break;
    }
    default:
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "The vector for the sorption coefficients must have 3 or 4 rows. Its a third degree polynomial; "
                "in case of 4 coefficients the constant term will be deleted");
        break;
    }
    }
}


void DataMoistureTransport::SetLastRelHumValue(double rLastRelHumValue)
{
    mLastRelHumValue = rLastRelHumValue;
}


void DataMoistureTransport::SetLastSorptionCoeff(FullVector<double,Eigen::Dynamic> rLastSorptionCoeff)
{
    switch (rLastSorptionCoeff.GetNumRows())
    {
    case 3:
    {
        mLastSorptionCoeff = rLastSorptionCoeff;
        break;
    }
    case 4:
    {
        mLastSorptionCoeff.Resize(3);
        for(int i=0; i<3; i++)
        {
            mLastSorptionCoeff(i) = rLastSorptionCoeff(i+1);
        }
        break;
    }
    default:
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__,
                "The vector for the sorption coefficients must have 3 or 4 rows. Its a third degree polynomial; "
                "in case of 4 coefficients the constant term will be deleted");
        break;
    }
    }
}


void DataMoistureTransport::SetDesorption(bool rSorptionHistoryDesorption)
{
    mSorptionHistoryDesorption = rSorptionHistoryDesorption;
}


bool DataMoistureTransport::operator==(const DataMoistureTransport& rhs) const
{
    bool value = this->mLastSorptionCoeff == rhs.mLastSorptionCoeff;
    value = value and (this->mLastRelHumValue == rhs.mLastRelHumValue);
    value = value and (this->mLastJunctionPoint == rhs.mLastJunctionPoint);
    value = value and (this->mCurrentJunctionPoint == rhs.mCurrentJunctionPoint);
    value = value and (this->mCurrentSorptionCoeff == rhs.mCurrentSorptionCoeff);
    value = value and (this->mLastSorptionCoeff == rhs.mLastSorptionCoeff);
    return value;
}


bool DataMoistureTransport::operator!=(const DataMoistureTransport& rhs) const
{
    return !this->operator==(rhs); 
}


double DataMoistureTransport::GetLastJunctionPoint() const
{
    return mLastJunctionPoint;
}


void DataMoistureTransport::SetLastJunctionPoint(double newLastJunctionPoint)
{
    mLastJunctionPoint = newLastJunctionPoint;
}


double DataMoistureTransport::GetCurrentJunctionPoint() const
{
    return mCurrentJunctionPoint;
}


void DataMoistureTransport::SetCurrentJunctionPoint(double newCurrentJunctionPoint)
{
    mCurrentJunctionPoint = newCurrentJunctionPoint;
}

namespace NuTo{namespace Constitutive{ namespace StaticData{
std::ostream& operator<<(std::ostream& os, const DataMoistureTransport& data)
{
    os << "SorptionHistoryDesorption: "     << data.mSorptionHistoryDesorption << "\n";
    os << "Last relative humidity value: "  << data.mLastRelHumValue << "\n";
    os << "Last junction point: "           << data.mLastJunctionPoint << "\n";
    os << "Current junction point: "        << data.mCurrentJunctionPoint << "\n";

    os << "Current sorption coefficients: " << data.mCurrentSorptionCoeff << "\n";
    os << "Last sorption coefficients: "    << data.mLastSorptionCoeff << "\n";
    return os;
}
}}} // namespaces

template <typename TStream>
void DataMoistureTransport::SerializeDataMoistureTransport(TStream &rStream)
{
    rStream.Serialize(mSorptionHistoryDesorption);
    rStream.Serialize(mLastRelHumValue);
    rStream.Serialize(mLastJunctionPoint);
    rStream.Serialize(mCurrentJunctionPoint);
    rStream.Serialize(mCurrentSorptionCoeff);
    rStream.Serialize(mLastRelHumValue);
}

template void DataMoistureTransport::SerializeDataMoistureTransport<NuTo::SerializeStreamIn>(SerializeStreamIn& rStream);
template void DataMoistureTransport::SerializeDataMoistureTransport<NuTo::SerializeStreamOut>(SerializeStreamOut& rStream);
