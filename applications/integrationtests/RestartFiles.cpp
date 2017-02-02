//
// Created by Thomas Titscher on 10/24/16.
//
#define BOOST_TEST_MODULE RestartFilesTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/tools/MeshGenerator.h"
#include "mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "mechanics/elements/ElementBase.h"

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}


namespace RestartFilesTest
{

void SetDummyStaticData(NuTo::Structure& rS, double rFactor)
{
    // set some static data values
    int gElementsTotal = rS.GroupGetElementsTotal();
    auto elementIds = rS.GroupGetMemberIds(gElementsTotal);
    for (unsigned int i = 0; i < elementIds.size(); ++i)
    {
        auto& e = *rS.ElementGetElementPtr(elementIds[i]);
        for (int ip = 0; ip < e.GetNumIntegrationPoints(); ++ip)
        {
            e.GetIPData().GetIPConstitutiveLaw(ip).GetData<NuTo::GradientDamageEngineeringStress>().SetData(0.1*i + ip*rFactor);
        }
    }
    rS.GroupDelete(gElementsTotal);
}

//! @brief creates a test structure
//! @param rS structure
//! @param rDummyValues true: fills it with some values; false: fills it with zeros
void CreateTestStructure(NuTo::Structure& rS, bool rDummyValues)
{
    rS.SetShowTime(false);
    rS.SetNumTimeDerivatives(1);
    int section = rS.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    int claw = rS.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    int interpolationType = rS.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    rS.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    rS.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    rS.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    NuTo::MeshGenerator::MeshRectangularPlane(rS, section, claw, interpolationType, {{2,2}}, {{42., 61.74}});
    rS.ElementTotalConvertToInterpolationType(); //VIRTHAMTODO :) I would expect that this method is called within MeshGenerator... right?

    // add a constraint --> dependent dof vector K
    rS.ConstraintLinearSetDisplacementNode(0, Eigen::Vector2d::UnitX(), 0.1337);

    rS.NodeBuildGlobalDofs();

    int gElementsTotal = rS.GroupGetElementsTotal();
    rS.ElementGroupAllocateAdditionalStaticData(gElementsTotal, 2);

    if (not rDummyValues)
        return;

    // set some node values
    auto nodeIds = rS.GroupGetMemberIds(rS.GroupGetNodesTotal());
    for (unsigned int i = 0; i < nodeIds.size(); ++i)
    {
        rS.NodeSetDisplacements(nodeIds[i], 0, Eigen::Vector2d({i+0.42e-7, i+0.6174}));
        rS.NodeSetDisplacements(nodeIds[i], 1, Eigen::Vector2d({i+0.42e-7, i+0.6174}));
    }
    // set some static data values
    SetDummyStaticData(rS, 1./17);
    rS.ElementTotalShiftStaticDataToPast();

    SetDummyStaticData(rS, 1./19);
    rS.ElementTotalShiftStaticDataToPast();

    SetDummyStaticData(rS, 1./21);
}

void CheckStaticData(NuTo::Structure& rSA, NuTo::Structure& rSB)
{
    // set some static data values
    int gElementsTotal = rSA.GroupGetElementsTotal();
    auto elementIds = rSA.GroupGetMemberIds(gElementsTotal);
    for (unsigned int i = 0; i < elementIds.size(); ++i)
    {
        auto& eA = *rSA.ElementGetElementPtr(elementIds[i]);
        auto& eB = *rSB.ElementGetElementPtr(elementIds[i]);
        for (int ip = 0; ip < eA.GetNumIntegrationPoints(); ++ip)
        {
            auto& dataA = eA.GetIPData().GetIPConstitutiveLaw(ip).GetData<NuTo::GradientDamageEngineeringStress>();
            auto& dataB = eB.GetIPData().GetIPConstitutiveLaw(ip).GetData<NuTo::GradientDamageEngineeringStress>();
            for (unsigned iData = 0; iData < dataA.GetNumData(); ++iData)
                BOOST_CHECK_CLOSE(dataA.GetData(iData), dataB.GetData(iData), 1.e-10);
        }
    }
}


void CheckDofs(NuTo::Structure& rSA, NuTo::Structure& rSB)
{
    int numTimeDerivatives = rSA.GetNumTimeDerivatives();
    for (int i = 0; i < numTimeDerivatives; ++i)
    {
        auto dofsA = rSA.NodeExtractDofValues(i);
        auto dofsB = rSB.NodeExtractDofValues(i);

        BOOST_CHECK_CLOSE(dofsA.J.Export().norm(), dofsB.J.Export().norm(), 1.e-10);
        BOOST_CHECK_CLOSE(dofsA.K.Export().norm(), dofsB.K.Export().norm(), 1.e-10);
    }
}

void NuToSerializeStructure(const std::string rFile, bool rIsBinary, NuTo::Structure& rA, NuTo::Structure& rB)
{
    CreateTestStructure(rA, true);

    { // write structure a --> file
        NuTo::SerializeStreamOut out(rFile, rIsBinary);
        out << rA;
    }

    CreateTestStructure(rB, false);

    { // read structure file --> b
        NuTo::SerializeStreamIn in(rFile, rIsBinary);
        in >> rB;
    }
}


BOOST_AUTO_TEST_CASE(RestartFiles_NuToSerializeStructureDofsText)
{
    NuTo::Structure a(2);
    NuTo::Structure b(2);
    NuToSerializeStructure("NuToSerializeStructureDofsText.dat" , false, a, b);
    CheckDofs(a, b);
}

BOOST_AUTO_TEST_CASE(RestartFiles_NuToSerializeStructureDofsBinary)
{
    NuTo::Structure a(2);
    NuTo::Structure b(2);
    NuToSerializeStructure("NuToSerializeStructureDofsBinary.dat" , true, a, b);
    CheckDofs(a, b);
}

BOOST_AUTO_TEST_CASE(RestartFiles_NuToSerializeStructureDataText)
{
    NuTo::Structure a(2);
    NuTo::Structure b(2);
    NuToSerializeStructure("NuToSerializeStructureDataText.dat" , false, a, b);
    CheckStaticData(a, b);
}

BOOST_AUTO_TEST_CASE(RestartFiles_NuToSerializeStructureDatasBinary)
{
    NuTo::Structure a(2);
    NuTo::Structure b(2);
    NuToSerializeStructure("NuToSerializeStructureDataBinary.dat" , true, a, b);
    CheckStaticData(a, b);
}

} // namespace RestartFilesTest
