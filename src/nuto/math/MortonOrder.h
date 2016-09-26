// $Id$
#pragma once

#include <stdint.h>
#include <limits>
namespace NuTo
{
//! @author Andrea Kessler, ISM
//! @date June 2013
//! @brief ... class for sorting in Morton order (z-order as space filling curve)
//! @brief ... inspired by 	from http://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/

class MortonOrder
{
public:
	static uint32_t EncodeMorton3D(uint32_t x,uint32_t y,uint32_t z)
	{
		return((Part1By2(z)<<2)+(Part1By2(y)<<1)+Part1By2(x));
	}

	// "Insert" a 0 bit after each of the 16 low bits of x
	static uint32_t Part1By1(uint32_t x)
	{
		x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
		x = (x ^ (x <<  8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
		x = (x ^ (x <<  4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
		x = (x ^ (x <<  2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
		x = (x ^ (x <<  1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
		return x;
	}

	// "Insert" two 0 bits after each of the 10 low bits of x
	static uint32_t Part1By2(uint32_t x)
	{
		x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
		x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
		return x;
	}
	// Inverse of Part1By1 - "delete" all odd-indexed bits
	static uint32_t Compact1By1(uint32_t x)
	{
		x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
		x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
		x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
		x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
		x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
		return x;
	}

	// Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
	static uint32_t Compact1By2(uint32_t x)
	{
		x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
		x = (x ^ (x >>  2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x ^ (x >>  4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x ^ (x >>  8)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
		return x;
	}

	static uint32_t DecodeMorton2X(uint32_t code)
	{
		return Compact1By1(code >> 0);
	}

	static uint32_t DecodeMorton2Y(uint32_t code)
	{
		return Compact1By1(code >> 1);
	}

	static uint32_t DecodeMorton3X(uint32_t code)
	{
		return Compact1By2(code >> 0);
	}

	static uint32_t DecodeMorton3Y(uint32_t code)
	{
		return Compact1By2(code >> 1);
	}

	static uint32_t DecodeMorton3Z(uint32_t code)
	{
		return Compact1By2(code >> 2);
	}
	//! @brief ... calculates neighbor of key shifted by +x+y+z
	static uint32_t Neighbor3D(uint32_t key,uint32_t dx,uint32_t dy,uint32_t dz)
	{
		uint32_t x=DecodeMorton3X(key);
		uint32_t y=DecodeMorton3Y(key);
		uint32_t z=DecodeMorton3Z(key);
		return EncodeMorton3D(x+dx,y+dy,z+dz);
		//return ((Part1By2(z+dz)<<2)+(Part1By2(y+dy)<<1)+Part1By2(x+dx));
	}
	//! @brief ... calculates neighbor of key shifted by +x+y+z
	static uint32_t NegNeighbor3D(uint32_t key,int32_t dx,int32_t dy,int32_t dz)
	{
		int32_t x=(int32_t) DecodeMorton3X(key);
		int32_t y=(int32_t) DecodeMorton3Y(key);
		int32_t z=(int32_t) DecodeMorton3Z(key);
		if (x+dx<0 || y+dy<0 || z+dz<0)
			return std::numeric_limits<uint32_t>::max (); // return max value, to have one neighbor which does not exist
		else
			return EncodeMorton3D( x+dx,y+dy,z+dz);
		//return ((Part1By2(z+dz)<<2)+(Part1By2(y+dy)<<1)+Part1By2(x+dx));
	}

};
} // namespace NuTo
