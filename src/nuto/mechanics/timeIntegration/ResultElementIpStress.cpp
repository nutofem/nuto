/*
 * ResultDispNode.cpp

 *
 *  Created on: Dec 18, 2013
 *      Author: junger
 */
#include <boost/assign/ptr_map_inserter.hpp>
#include "nuto/mechanics/elements/ElementOutputIpData.h"

#include "nuto/mechanics/timeIntegration/ResultElementIpStress.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::ResultElementIpStress::ResultElementIpStress(const std::string& rIdent, int rElementId) : ResultElementIpBase(rIdent, rElementId)
{
}

//! @brief calculate the relevant nodal dofs
void NuTo::ResultElementIpStress::CalculateValues(const StructureBase& rStructure, NuTo::FullMatrix<double, 1, Eigen::Dynamic>& rValues)const
{
	const ElementBase*  element(rStructure.ElementGetElementPtr(mElementId));

    //determine the ipdata and determine the map
	boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase> elementOutput;
	boost::assign::ptr_map_insert<ElementOutputIpData>( elementOutput )( Element::IP_DATA ,IpData::ENGINEERING_STRESS);

    //calculate the element solution (this can't be const, since update is also a possible step in evaluate which is not constant)
	const_cast<ElementBase*>(element)->Evaluate(elementOutput);

    //assign the outputs
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>* engineeringStress(0);
    for (auto itElementOutput=elementOutput.begin(); itElementOutput!=elementOutput.end(); itElementOutput++)
    {
        switch (itElementOutput->second->GetIpDataType())
        {
		case NuTo::IpData::ENGINEERING_STRESS:
			engineeringStress = &(itElementOutput->second->GetFullMatrixDouble());
		break;
		default:
			throw MechanicsException("[NuTo::ElementBase::Visualize] other ipdatatypes not supported.");
        }
    }

    //iterate over all ips
    assert(engineeringStress->GetNumColumns()==element->GetNumIntegrationPoints());
    assert(engineeringStress->GetNumRows()==6);
    for (int count=0; count<engineeringStress->GetNumColumns(); count++)
    {
    	rValues.SetBlock(0,count*6,engineeringStress->GetColumn(count).Trans());
    }
    return;
}

//! @brief number of data points per time step (e.g. number of displacement components of a node
int NuTo::ResultElementIpStress::GetNumData(const StructureBase& rStructure)const
{
	const ElementBase* element(rStructure.ElementGetElementPtr(mElementId));
	return element->GetNumIntegrationPoints()*6;
}


