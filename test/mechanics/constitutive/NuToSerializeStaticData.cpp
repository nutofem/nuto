#define BOOST_TEST_MODULE NuToSerializeStaticDataTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <eigen3/Eigen/Core>

#include "nuto/base/serializeStream/SerializeStreamOut.h"
#include "nuto/base/serializeStream/SerializeStreamIn.h"

#include "nuto/mechanics/constitutive/staticData/Leaf.h"
#include "nuto/mechanics/constitutive/staticData/Composite.h"

// needed for building with clang when boost test has been built with gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size()));
}


namespace NuToSerializeStaticData
{

class A
{
public:

    A(const Eigen::Matrix2d& rData) : mData(rData) {}

    void NuToSerializeSave(NuTo::SerializeStreamOut& rStream)
    {
        rStream.Serialize(mData);
    }

    void NuToSerializeLoad(NuTo::SerializeStreamIn& rStream)
    {
        rStream.Serialize(mData);
    }

    Eigen::Matrix2d mData;

};

template <typename T>
T SerializeLeaf(const std::string& rFile, bool rIsBinary, T& rValue, T rZero)
{
    auto leaf = NuTo::Constitutive::StaticData::Leaf<T>::Create(rValue);
    {
        NuTo::SerializeStreamOut s(rFile, rIsBinary);
        s << *leaf;
    }

    auto leafFromFile = NuTo::Constitutive::StaticData::Leaf<T>::Create(rZero);

    {
        NuTo::SerializeStreamIn s(rFile, rIsBinary);
        s >> *leafFromFile;
    }
    T data = leafFromFile->GetData();
    delete leaf;
    delete leafFromFile;
    return data;
}

BOOST_AUTO_TEST_CASE(NuToSerializeLeafDouble)
{
    double d = 5.;
    BOOST_CHECK_CLOSE(d, SerializeLeaf("LeafDoubleText.dat", false, d, 0.), 1.e-10);
    BOOST_CHECK_CLOSE(d, SerializeLeaf("LeafDoubleBinary.dat", true, d, 0.), 1.e-10);
}


template <typename T>
T SerializeLeafBase(const std::string& rFile, bool rIsBinary, T& rValue, T rZero)
{
    NuTo::Constitutive::StaticData::Component* base = NuTo::Constitutive::StaticData::Leaf<T>::Create(rValue);
    {
        NuTo::SerializeStreamOut s(rFile, rIsBinary);
        s << *base;
    }

    NuTo::Constitutive::StaticData::Component* baseFromFile = NuTo::Constitutive::StaticData::Leaf<T>::Create(rZero);

    {
        NuTo::SerializeStreamIn s(rFile, rIsBinary);
        s >> *baseFromFile;
    }
    T data = dynamic_cast<NuTo::Constitutive::StaticData::Leaf<T>&>(*baseFromFile).GetData();
    delete base;
    delete baseFromFile;
    return data;
}


BOOST_AUTO_TEST_CASE(NuToSerializeLeafDoubleBase)
{
    double d = 5.;
    BOOST_CHECK_CLOSE(d, SerializeLeafBase("LeafBaseDoubleText.dat", false, d, 0.), 1.e-10);
    BOOST_CHECK_CLOSE(d, SerializeLeafBase("LeafBaseDoubleBinary.dat", true, d, 0.), 1.e-10);
}

BOOST_AUTO_TEST_CASE(NuToSerializeLeafABase)
{
    A a(Eigen::Matrix2d::Random());
    A zero(Eigen::Matrix2d::Zero(2,2));
    BOOST_CHECK_CLOSE(a.mData.norm(), SerializeLeafBase("LeafBaseAText.dat", false, a, zero).mData.norm(), 1.e-10);
    BOOST_CHECK_CLOSE(a.mData.norm(), SerializeLeafBase("LeafBaseABinary.dat", true, a, zero).mData.norm(), 1.e-10);
}


NuTo::Constitutive::StaticData::Component* GetComposite(double rA, double rB, double rC)
{
    using namespace NuTo::Constitutive::StaticData;
    auto composite = Composite::Create();
    composite->AddComponent(Leaf<double>::Create(rA));
    composite->AddComponent(Composite::Create());
    auto& subComposite = dynamic_cast<Composite&>(composite->GetComponent(1));
    subComposite.AddComponent(Leaf<double>::Create(rB));
    subComposite.AddComponent(Leaf<double>::Create(rC));

    composite->AllocateAdditionalData(3);

    return composite;
}

void CheckLeafDouble(NuTo::Constitutive::StaticData::Component& c1, NuTo::Constitutive::StaticData::Component& c2)
{
    auto& l1 = dynamic_cast<NuTo::Constitutive::StaticData::Leaf<double>&>(c1);
    auto& l2 = dynamic_cast<NuTo::Constitutive::StaticData::Leaf<double>&>(c2);

    BOOST_CHECK_EQUAL(l1.GetNumData(), l2.GetNumData());

    for (auto i = 0; i < l1.GetNumData(); ++i)
        BOOST_CHECK_CLOSE(l1.GetData(i), l2.GetData(i), 1.e-10);
}

void CheckComposite(NuTo::Constitutive::StaticData::Component& c1, NuTo::Constitutive::StaticData::Component& c2)
{
    using namespace NuTo::Constitutive::StaticData;
    auto& c1composite = dynamic_cast<Composite&>(c1);
    auto& c2composite = dynamic_cast<Composite&>(c2);

    CheckLeafDouble(c1composite.GetComponent(0), c2composite.GetComponent(0));

    auto& c1subcomposite = dynamic_cast<Composite&>(c1composite.GetComponent(1));
    auto& c2subcomposite = dynamic_cast<Composite&>(c2composite.GetComponent(1));

    CheckLeafDouble(c1subcomposite.GetComponent(0), c2subcomposite.GetComponent(0));
    CheckLeafDouble(c1subcomposite.GetComponent(1), c2subcomposite.GetComponent(1));
}

void SerializeComposite(const std::string& rFile, bool rIsBinary)
{
    auto composite = GetComposite(1., 2., 3.);
    auto compositeFromFile = GetComposite(0., 0., 0.);

    {
        NuTo::SerializeStreamOut s(rFile, rIsBinary);
        s << *composite;
    }

    {
        NuTo::SerializeStreamIn s(rFile, rIsBinary);
        s >> *compositeFromFile;
    }
    CheckComposite(*composite, *compositeFromFile);
    delete composite;
    delete compositeFromFile;
}



BOOST_AUTO_TEST_CASE(SerializeCompositeBase)
{
    SerializeComposite("ComponentText", false);
    SerializeComposite("ComponentBinary", true);
}

} // namespace




