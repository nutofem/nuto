/*
 * InputReader.h
 *
 *  Created on: 17 Mar 2014
 *      Author: ttitsche
 */

#ifndef INPUTREADER_H_
#define INPUTREADER_H_

#include "nuto/math/FullMatrix.h"
#include "nuto/math/FullVector.h"


#include <string>
#include <fstream>


/* Example File:

particle simulation input file

 ------------------------------------------
| simulation parameters and stop criterion |
 ------------------------------------------
numEventsMax = 1.e10;
timeMax = 10.;
timePrintOut = 1.;
initialTimeBarrier = 1.;
randomVelocityRange = 0.;
growthRates = .1;
numThreads = 4;


 ------------------------------------------
|        bounding box definition           |
 ------------------------------------------
boxType = 0;
boundingBox = 0., 100.,
              0., 100.,
              0., 100.;

 ------------------------------------------
|           particle definition            |
 ------------------------------------------
definitionByGradingCurve = true;

numSieves = 3;
   sieve1 =  8., 16., 0.40;
   sieve2 =  4.,  8., 0.24;
   sieve3 =  2.,  4., 0.15;
volumeFraction = 0.70;
absoluteDistance = 0.0;

 */


namespace NuTo
{

class InputReader
{
public:
	InputReader(std::string rFileName);
	void PrintInput();
	double GetAbsoluteDistance() const;
	const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetBoundingBox() const;
	int GetBoxType() const;
	const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetGradingCurve() const;
	double GetGrowthRates() const;
	double GetInitialTimeBarrier() const;
	bool isIsDefinedByGradingCurve() const;
	long GetNumEventsMax() const;
	int GetNumParticles() const;
	double GetRandomVelocityRange() const;
	double GetTimeMax() const;
	double GetTimePrintOut() const;
	double GetVolumeFraction() const;
	int GetNumThreads() const;
	const std::string& GetDirectory() const;

	double GetShrinkage() const
	{
		return mShrinkage;
	}

private:

	void OpenFile(std::string rFileName);
	void SkipToNextData();
	double ReadNumber();
	FullVector<double, Eigen::Dynamic> ReadVector();
	bool ReadBool();
	std::string ReadString();

	void ReadFile();
	void ReadSimulationParameters();
	void ReadBoundingBox();
	void ReadGradingCurve();
	void ReadMonodisperse();

	std::ifstream mFile;

	std::string mDirectory;

	long mNumEventsMax;
	double mTimeMax;
	double mTimePrintOut;
	double mInitialTimeBarrier;
	double mRandomVelocityRange;
	double mGrowthRates;
	int mNumThreads;

	int mBoxType;
	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mBoundingBox;

	bool mIsDefinedByGradingCurve;

	FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mGradingCurve;
	double mVolumeFraction;
	double mAbsoluteDistance;

	int mNumParticles;

	double mShrinkage;

};

} /* namespace NuTo */
#endif /* INPUTREADER_H_ */
