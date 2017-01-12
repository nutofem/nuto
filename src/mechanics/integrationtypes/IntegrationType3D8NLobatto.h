// $Id$
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType3D8NLobatto_Def.h"

namespace NuTo
{
    //! @brief constructor
    template <int T>
	IntegrationType3D8NLobatto<T>::IntegrationType3D8NLobatto()
	{
	    mCoordinates.resize(T*T*T);
	    mWeights.resize(T*T*T);

	    NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;
	    NuTo::IntegrationType1D2NLobatto4Ip Lobatto1D2N4Ip;
	    NuTo::IntegrationType1D2NLobatto5Ip Lobatto1D2N5Ip;

	    NuTo::IntegrationType1D* integrationType1D(0);
		switch(T)
		{
	    case 3:
	    {
	    	integrationType1D = &Lobatto1D2N3Ip;
	    	break;
	    }
	    case 4:
	    {
	    	integrationType1D = &Lobatto1D2N4Ip;
	    	break;
	    }
	    case 5:
	    {
	    	integrationType1D = &Lobatto1D2N5Ip;
	    	break;
	    }
	    default:
	    	throw MechanicsException("[IntegrationType3D8NLobatto<T>::IntegrationType3D8NLobatto] Only implemented for orders 3,4 and 5");
		}

		assert (integrationType1D->GetNumIntegrationPoints()==T);

	    std::array<double,T> coordinates1D;
	    std::array<double,T> weights1D;

	    // get the 1D integration point coordinates and weights
	    for (int i = 0; i < T; i++)
	    {
			coordinates1D[i] = integrationType1D->GetLocalIntegrationPointCoordinates(i)[0];
	        weights1D[i] = integrationType1D->GetIntegrationPointWeight(i);
	    }

	    // calculate the 2D integratration point coordinates and weights
	    int ipNum = 0;
	    for (int iz = 0; iz < T; iz++)
	        for (int iy= 0; iy < T; iy++)
		        for (int ix= 0; ix < T; ix++, ipNum++)
	        {
	            mWeights[ipNum] = weights1D[ix]*weights1D[iy]*weights1D[iz];
	            mCoordinates[ipNum][0] = coordinates1D[ix];
	            mCoordinates[ipNum][1] = coordinates1D[iy];
	            mCoordinates[ipNum][2] = coordinates1D[iz];
//	            std::cout << "ip " << ipNum << " with coord " <<  mCoordinates[ipNum][0] << " " << mCoordinates[ipNum][1] << " " << mCoordinates[ipNum][2] << " and weight " << mWeights[ipNum] << std::endl;
	        }
	}

#ifndef SWIG
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <int T>
    template<class Archive>
    void IntegrationType3D8NLobatto<T>::serialize(Archive & ar, const unsigned int version)
    {}
#endif // ENABLE_SERIALIZATION
#endif //SWIG

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    template <int T>
    Eigen::VectorXd IntegrationType3D8NLobatto<T>::GetLocalIntegrationPointCoordinates(int rIpNum) const
	{
        assert(rIpNum>=0 && rIpNum<T*T*T);
        return mCoordinates[rIpNum];
	}


    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    template <int T>
    int IntegrationType3D8NLobatto<T>::GetNumIntegrationPoints()const
	{
    	return T*T*T;
	}

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    template <int T>
    double IntegrationType3D8NLobatto<T>:: GetIntegrationPointWeight(int rIpNum)const
	{
        assert(rIpNum>=0 && rIpNum<T*T*T);
    	return mWeights[rIpNum];
	}

    //! @brief returns a string with the identifier of the integration type
    //! @return identifier
    template <int T>
	eIntegrationType IntegrationType3D8NLobatto<T>::GetEnumType() const
    {
    	switch(T)
		{
	    case 3: return eIntegrationType::IntegrationType3D8NLobatto3x3x3Ip;
	    case 4: return eIntegrationType::IntegrationType3D8NLobatto4x4x4Ip;
	    case 5: return eIntegrationType::IntegrationType3D8NLobatto5x5x5Ip;
	    default:
	    	throw MechanicsException("[IntegrationType3D8NLobatto<T>::IntegrationType3D8NLobatto] Only implemented for orders 3,4 and 5");
		}
    }

#ifdef ENABLE_VISUALIZE
    template <int T>
    void IntegrationType3D8NLobatto<T>::GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const
	{


    	NumVisualizationPoints = (T+1)*(T+1)*(T+1);

    	double x,y,z;
        for (int countZ = 0; countZ<=T;countZ++)
        {
        	switch (countZ)
        	{
        	case 0:
        		z=-1;
        		break;
        	case T:
        		z=1.;
        		break;
        	default:
        		z=0.5*(mCoordinates[countZ][0]+mCoordinates[countZ-1][0]);
        	}
            for (int countY = 0; countY<=T;countY++)
            {
            	switch (countY)
            	{
            	case 0:
            		y=-1;
            		break;
            	case T:
            		y=1.;
            		break;
            	default:
            		y=0.5*(mCoordinates[countY][0]+mCoordinates[countY-1][0]);
            	}
                for (int countX = 0; countX<=T;countX++)
                {
                	switch (countX)
                	{
                	case 0:
                		x=-1;
                		break;
                	case T:
                		x=1.;
                		break;
                	default:
                		x=0.5*(mCoordinates[countX][0]+mCoordinates[countX-1][0]);
                	}

					VisualizationPointLocalCoordinates.push_back(x);
					VisualizationPointLocalCoordinates.push_back(y);
					VisualizationPointLocalCoordinates.push_back(z);
                }
            }
        }

        NumVisualizationCells = T*T*T;
        int theIp(0);
        for (int countZ = 0; countZ<T;countZ++)
        {
            for (int countY = 0; countY<T;countY++)
            {
                for (int countX = 0; countX<T;countX++,theIp++)
                {
                    VisualizationCellType.push_back(NuTo::eCellTypes::HEXAHEDRON);
                    VisualizationCellsIncidence.push_back(countZ*(T+1)*(T+1)+countY*(T+1)+countX);
                    VisualizationCellsIncidence.push_back(countZ*(T+1)*(T+1)+countY*(T+1)+countX+1);
                    VisualizationCellsIncidence.push_back(countZ*(T+1)*(T+1)+(countY+1)*(T+1)+countX+1);
                    VisualizationCellsIncidence.push_back(countZ*(T+1)*(T+1)+(countY+1)*(T+1)+countX);
                    VisualizationCellsIncidence.push_back((countZ+1)*(T+1)*(T+1)+countY*(T+1)+countX);
                    VisualizationCellsIncidence.push_back((countZ+1)*(T+1)*(T+1)+countY*(T+1)+countX+1);
                    VisualizationCellsIncidence.push_back((countZ+1)*(T+1)*(T+1)+(countY+1)*(T+1)+countX+1);
                    VisualizationCellsIncidence.push_back((countZ+1)*(T+1)*(T+1)+(countY+1)*(T+1)+countX);
                    VisualizationCellsIP.push_back(theIp);
                }
            }
        }
	}
#endif // ENABLE_VISUALIZE

}

