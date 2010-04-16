// $Id$

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

// create a new section
void NuTo::StructureBase::SectionCreate(const std::string& rIdent, const std::string& rType)
{
    // convert section type string to upper case
    std::string SectionTypeString;
    std::transform(rType.begin(), rType.end(), std::back_inserter(SectionTypeString), (int(*)(int)) toupper);

    // get section type from string
    SectionBase::eSectionType SectionType;
    if (SectionTypeString == "TRUSS")
    {
        SectionType = SectionBase::TRUSS;
    }
    else if (SectionTypeString == "PLANE_STRAIN")
    {
        SectionType = SectionBase::PLANE_STRAIN;
    }
    else if (SectionTypeString == "PLANE_STRESS")
    {
        SectionType = SectionBase::PLANE_STRESS;
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionCreate] invalid section type.");
    }
    this->SectionCreate(rIdent, SectionType);
}

// create a new section
void NuTo::StructureBase::SectionCreate(const std::string& rIdent, SectionBase::eSectionType rType)
{
    // check if section identifier exists
    boost::ptr_map<std::string,SectionBase>::iterator it = this->mSectionMap.find(rIdent);
    if (it == this->mSectionMap.end())
    {
        // create new section
        SectionBase* SectionPtr;
        switch (rType)
        {
        case SectionBase::TRUSS:
            SectionPtr = new SectionTruss();
            break;
        case SectionBase::PLANE_STRAIN:
        case SectionBase::PLANE_STRESS:
            SectionPtr = new SectionPlane(rType);
            break;
        default:
            throw NuTo::MechanicsException("[NuTo::StructureBase::SectionCreate] invalid section type.");
        }

        // add section to map (insert does not allow const keys!!!!)
        this->mSectionMap.insert(const_cast<std::string&>(rIdent), SectionPtr);
    }
    else
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionCreate] Section already exists.");
    }
}

// delete an existing section
void NuTo::StructureBase::SectionDelete(const std::string& rIdent)
{
    // find section identifier in map
    boost::ptr_map<std::string,SectionBase>::iterator it = this->mSectionMap.find(rIdent);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionDelete] Section does not exist.");
    }
    else
    {
        this->mSectionMap.erase(it);
    }
}

// get section pointer from section identifier
NuTo::SectionBase* NuTo::StructureBase::SectionGetSectionPtr(const std::string& rIdent)
{
    boost::ptr_map<std::string,SectionBase>::iterator it = this->mSectionMap.find(rIdent);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionGetSectionPtr] Section does not exist.");
    }
    return it->second;
}

// get section pointer from section identifier
const NuTo::SectionBase* NuTo::StructureBase::SectionGetSectionPtr(const std::string& rIdent) const
{
    boost::ptr_map<std::string,SectionBase>::const_iterator it = this->mSectionMap.find(rIdent);
    if (it == this->mSectionMap.end())
    {
        throw NuTo::MechanicsException("[NuTo::StructureBase::SectionGetSectionPtr] Section does not exist.");
    }
    return it->second;
}

// get section identifier from section pointer
std::string NuTo::StructureBase::SectionGetId(const NuTo::SectionBase* rSectionPtr) const
{
    for (boost::ptr_map<std::string,SectionBase>::const_iterator it = mSectionMap.begin(); it!= mSectionMap.end(); it++)
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
    for (boost::ptr_map<std::string,SectionBase>::const_iterator it = mSectionMap.begin(); it!= mSectionMap.end(); it++)
    {
        std::cout << "  Section: " << it->first << std::endl;
        it->second->Info(rVerboseLevel);
    }
}

void NuTo::StructureBase::SectionInfo(const std::string& rIdent, unsigned short rVerboseLevel) const
{
    const SectionBase* SectionPtr = this->SectionGetSectionPtr(rIdent);
    std::cout << "  Section: " << rIdent << std::endl;
    SectionPtr->Info(rVerboseLevel);
}

// set section cross-section area
void NuTo::StructureBase::SectionSetArea(const std::string& rIdent, double rArea)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rIdent);
        SectionPtr->SetArea(rArea);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetArea] error setting section cross-section area.");
        throw e;
    }
}

// get section cross-section area
double NuTo::StructureBase::SectionGetArea(const std::string& rIdent) const
{
    double area;
    try
    {
        const SectionBase* SectionPtr = this->SectionGetSectionPtr(rIdent);
        area = SectionPtr->GetArea();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionGetArea] error getting section cross-section area.");
        throw e;
    }
    return area;
}

// set section thickness
void NuTo::StructureBase::SectionSetThickness(const std::string& rIdent, double rThickness)
{
    try
    {
        SectionBase* SectionPtr = this->SectionGetSectionPtr(rIdent);
        SectionPtr->SetThickness(rThickness);
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionSetThickness] error setting section thickness.");
        throw e;
    }
}

// get section thickness
double NuTo::StructureBase::SectionGetThickness(const std::string& rIdent) const
{
    double thickness;
    try
    {
        const SectionBase* SectionPtr = this->SectionGetSectionPtr(rIdent);
        thickness = SectionPtr->GetThickness();
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureBase::SectionGetThickness] error getting section thickness.");
        throw e;
    }
    return thickness;
}
