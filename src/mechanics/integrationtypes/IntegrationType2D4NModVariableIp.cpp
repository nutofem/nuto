// $Id$

#include <assert.h>

#include <cmath>
#include <sstream>
#include <iostream>

#include <boost/foreach.hpp>

#include "mechanics/integrationtypes/IntegrationType2D4NModVariableIp.h"
#include "mechanics/integrationtypes/IntegrationType2D4NConstVariableIp.h"
#include "mechanics/integrationtypes/IntegrationPointBase.h"
#include "mechanics/MechanicsException.h"

//! @brief constructor
//! @author Daniel Arnold, ISM
//! @date February 2011
NuTo::IntegrationType2D4NModVariableIp::IntegrationType2D4NModVariableIp(const std::string rName,const int rNumIp)
{
	mName=rName;
	const NuTo::IntegrationType2D4NConstVariableIp orgIPScheme(rNumIp);

	double locCoords[2]={0,0};
	double weight=0.0;
	int numIp=orgIPScheme.GetNumIntegrationPoints();
	const double width=2./sqrt((double)numIp);
	const double height=2./sqrt((double)numIp);
	for(int ip=0 ; ip<numIp;++ip)
	{
		//! get the original weighting factor
		weight=orgIPScheme.GetIntegrationPointWeight(ip);
		//! get the original local coordinates of this IP
		orgIPScheme.GetLocalIntegrationPointCoordinates2D(ip,locCoords);
		std::vector< double > locCoordsVec(0);
		BOOST_FOREACH(double d, locCoords)
			locCoordsVec.push_back(d);
		//! build up the bounding box
		std::vector< double > thisBoundingBox(8);
			thisBoundingBox[0]=locCoords[0]-0.5*width;
			thisBoundingBox[1]=locCoords[1]-0.5*height;
			thisBoundingBox[2]=locCoords[0]+0.5*width;
			thisBoundingBox[3]=locCoords[1]-0.5*height;
			thisBoundingBox[4]=locCoords[0]+0.5*width;
			thisBoundingBox[5]=locCoords[1]+0.5*height;
			thisBoundingBox[6]=locCoords[0]-0.5*width;
			thisBoundingBox[7]=locCoords[1]+0.5*height;
		//! with this informations we build up the new IP
		IntegrationPointBase* ptrIP(0);
		ptrIP = new IntegrationPointBase(locCoordsVec, weight, thisBoundingBox);
		//! and finally insert it into the map
		mIpMap.insert(ip,ptrIP);
	}
}

