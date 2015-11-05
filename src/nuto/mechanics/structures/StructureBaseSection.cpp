// $Id$

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionSpring.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/sections/SectionVolume.h"
#include "nuto/mechanics/sections/SectionFibreMatrixBond.h"

// create a new section
int NuTo::StructureBase::SectionCreate(const std::string& rType)
{
    // convert section type string to upper case
    std::string SectionTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(SectionTypeString), (int (*)(int)) toupper);

    // get section type from string
Section    ::eSectionType SectionType;
    if (SectionTypeString == "TRUSS")
    {
        SectionType = Section::TRUSS;
    } else if (SectionTypeString == "PLANE_STRAIN")
    {
        SectionType = Section::PLANE_STRAIN;
    } else if (SectionTypeString == "PLANE_STRESS")
    {
        SectionType = Section::PLANE_STRESS;
    } else if (SectionTypeString == "VOLUME")
    {
        SectionType = Section::VOLUME;
    } else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionCreate] invalid section type.");
    }
    return this->SectionCreate(SectionType);
}

// create a new section
int NuTo::StructureBase::SectionCreate(Section::eSectionType rType)
{
    //find unused integer id
    int id(0);
    boost::ptr_map<int, SectionBase>::iterator it = mSectionMap.find(id);
    while (it != mSectionMap.end())
    {
        id++;
        it = mSectionMap.find(id);
    }

    // create new section
    SectionBase* SectionPtr;
    switch (rType)
    {
    case Section::SPRING:
        SectionPtr = new SectionSpring();
        break;
    case Section::TRUSS:
        SectionPtr = new SectionTruss();
        break;
    case Section::PLANE_STRAIN:
    case Section::PLANE_STRESS:
        SectionPtr = new SectionPlane(rType);
        break;
    case Section::VOLUME:
        SectionPtr = new SectionVolume();
        break;
    case Section::FIBRE_MATRIX_BOND:
        SectionPtr = new SectionFibreMatrixBond();
        break;
    default:
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionCreate] invalid section type.");
    }

    // add section to map (insert does not allow const keys!!!!)
    this->mSectionMap.insert(id, SectionPtr);
    return id;
}

// delete an existing section
void NuTo::StructureBase::SectionDelete(int rId)
{
    // find section identifier in map
    boost::ptr_map<int, SectionBase>::iterator it = this->mSectionMap.find(rId);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionDelete] Section does not exist.");
    } else
    {
        this->mSectionMap.erase(it);
    }
}

// get section pointer from section identifier
NuTo::SectionBase* NuTo::StructureBase::SectionGetSectionPtr(int rId)
{
    boost::ptr_map<int, SectionBase>::iterator it = this->mSectionMap.find(rId);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionGetSectionPtr] Section does not exist.");
    }
    return it->second;
}

// get section pointer from section identifier
const NuTo::SectionBase* NuTo::StructureBase::SectionGetSectionPtr(int rId) const
{
    boost::ptr_map<int, SectionBase>::const_iterator it = this->mSectionMap.find(rId);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionGetSectionPtr] Section does not exist.");
    }
    return it->second;
}

// get section identifier from section pointer
int NuTo::StructureBase::SectionGetId(const NuTo::SectionBase* rSectionPtr) const
{
    for (boost::ptr_map<int, SectionBase>::const_iterator it = mSectionMap.begin(); it != mSectionMap.end(); it++)
    {
        if (it->second == rSectionPtr)
        {
            return it->first;
        }
    }
    throw MechanicsException("[NuTo::StructureBase::GetSectionId] Section does not exist.");
}

// info routines
void NuTo::StructureBase::SectionInfo(unsigned short rVerboseLevel) const
{
    std::cout << "Number of sections: " << this->GetNumSections() << std::endl;
    for (boost::ptr_map<int, SectionBase>::const_iterator it = mSectionMap.begin(); it != mSectionMap.end(); it++)
    {
        std::cout << "  Section: " << it->first << std::endl;
        it->second->Info(rVerboseLevel);
    }
}

void NuTo::StructureBase::SectionInfo(int rId, unsigned short rVerboseLevel) const
{
    const SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
    std::cout << "  Section: " << rId << std::endl;
    SectionPtr->Info(rVerboseLevel);
}

// set section cross-section area
void NuTo::StructureBase::SectionSetArea(int rId, double rArea)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        SectionPtr->SetArea(rArea);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetArea] error setting section cross-section area.");
        throw e;
    }
}

// get section cross-section area
double NuTo::StructureBase::SectionGetArea(int rId) const
{
    double area;
    try
    {
        const SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        area = SectionPtr->GetArea();
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionGetArea] error getting section cross-section area.");
        throw e;
    }
    return area;
}

// set section thickness
void NuTo::StructureBase::SectionSetThickness(int rId, double rThickness)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        SectionPtr->SetThickness(rThickness);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetThickness] error setting section thickness.");
        throw e;
    }
}

// get section thickness
double NuTo::StructureBase::SectionGetThickness(int rId) const
{
    double thickness;
    try
    {
        const SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        thickness = SectionPtr->GetThickness();
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionGetThickness] error getting section thickness.");
        throw e;
    }
    return thickness;
}

void NuTo::StructureBase::SectionSetCircumference(int rId, double rCircumference)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        SectionPtr->SetCircumference(rCircumference);
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(std::string(__PRETTY_FUNCTION__) + ":\t error setting section circumference.");
        throw e;
    }
}

double NuTo::StructureBase::SectionGetCircumference(int rId) const
{
    double circumference;
    try
    {
        const SectionBase* SectionPtr = this->SectionGetSectionPtr(rId);
        circumference = SectionPtr->GetCircumference();
    } catch (NuTo::MechanicsException& e)
    {
        e.AddMessage(std::string(__PRETTY_FUNCTION__) + ":\t error getting section circumference.");
        throw e;
    }
    return circumference;
}

