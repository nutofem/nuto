// $Id$
#ifndef NODE_KEY_H
#define NODE_KEY_H

#include "nuto/mechanics/nodes/NodeEnum.h"

namespace NuTo
{
class NodeKey
{
	public:
	short mNumCoord;
	short mNumTimeDerivative;
	short mNumDisp;
	short mNumFineScaleDisp;
	short mNumRotations;
	short mNumTemperature;
	bool  mNonlocalData;
	bool  mRadius;
	bool  mGrid;

	NodeKey(short rNumCoord, short rNumTimeDerivative, short rNumDisp, short rNumFineScaleDisp=0, short rNumRotations=0, short rNumTemperature=0,
			bool rNonlocalData=false, bool rRadius=false, bool rGrid=false)
	: mNumCoord(rNumCoord), mNumTimeDerivative(rNumTimeDerivative), mNumDisp(rNumDisp), mNumFineScaleDisp(rNumFineScaleDisp),
	  mNumRotations(rNumRotations), mNumTemperature(rNumTemperature), mNonlocalData(rNonlocalData), mRadius(rRadius), mGrid(rGrid)
	{
	}

	bool operator<(const NodeKey &right) const
	{
		if ( mNumCoord == right.mNumCoord )
		{
			if ( mNumTimeDerivative == right.mNumTimeDerivative )
			{
				if ( mNumDisp == right.mNumDisp )
				{
					if ( mNumFineScaleDisp == right.mNumFineScaleDisp )
					{
						if ( mNumRotations == right.mNumRotations )
						{
							if ( mNumTemperature == right.mNumTemperature )
							{
								if ( mNonlocalData == right.mNonlocalData )
								{
									if ( mRadius == right.mRadius )
									{
										return mGrid < right.mGrid;
									}
									else
									{
										return mRadius < right.mNonlocalData;
									}
								}
								else
								{
									return mNonlocalData < right.mNonlocalData;
								}
							}
							else
							{
								return mNumTemperature < right.mNumTemperature;
							}
						}
						else
						{
							return mNumRotations < right.mNumRotations;
						}
					}
					else
					{
						return mNumFineScaleDisp < right.mNumFineScaleDisp;
					}
				}
				else
				{
					return mNumDisp < right.mNumDisp;
				}
			}
			else
			{
				return mNumTimeDerivative < right.mNumTimeDerivative;
			}
		}
		else
		{
			return mNumCoord < right.mNumCoord;
		}
	}
};

class MapNodeKeyNodeEnum
{
public:
	MapNodeKeyNodeEnum();
	~MapNodeKeyNodeEnum();
	NuTo::Node::eNodeType GetEnum(const NodeKey& rKey)const;

private:
	static std::map<NodeKey, NuTo::Node::eNodeType> mMapNodeKeyNodeEnum;
};
}//namespace NuTo

#endif //NODE_KEY_H
