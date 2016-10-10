/*
 * InputReader.h
 *
 *  Created on: 17 Mar 2014
 *      Author: ttitsche
 */

#pragma once

#include "nuto/math/FullMatrix_Def.h"


#include <string>
#include <fstream>

/* Example File:

particle simulation input file

1) avoid equal signes in comments
2) avoid comments before ';' in vector definitions
 ------------------------------------------
| simulation parameters and stop criterion |
 ------------------------------------------
directory = myNutoExamples/c++/GradingCurveA/;

numEventsMax = 1.e10;              Simulation stops after [numEventsMax] events
timeMax = 600.;                    Simulation stops after [timeMax] seconds
timePrintOut = 2.;                 There is a print out every [timePrintOut] seconds
initialTimeBarrier = .1;           Adaption during simulation, not that important
randomVelocityRange = 0.1;         equal distr. range[randomVelocityRange],center[0.0]
growthRates = 1.;                  Growth rate for all particles
numThreads = 8;                    Number of threads for parallelization

 ------------------------------------------
|        bounding box definition           |
 ------------------------------------------
boxType = 0;                       0.. cube | 1.. dogbone | 2.. cylinder
boundingBox = 0.,80.,
              0.,80.,
              0.,80.;            x-y-z- range
is2D = false;
 ------------------------------------------
|           particle definition            |
 ------------------------------------------

numSieves = 3;                     size of the grading curve
   sieve1 =  8., 16., 0.40;        1st row: min particle size
   sieve2 =  4.,  8., 0.24;        2nd row: max particle size
   sieve3 =  2.,  4., 0.15;        3rd row: volume fraction of this class


volumeFraction = 0.7;             overall particle volume fraction
absoluteDistance = 0.0;            minimal distance between two particles

shrinkage = 0.10;              0.0 -> no EDMD

*/


namespace NuTo
{
template <class T, int rows> class FullVector;
class InputReader
{
public:
	InputReader(std::string rFileName);
	void PrintInput();
	double GetAbsoluteDistance() const;
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetBoundingBox() const;
	int GetTypeOfSpecimen() const;
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetGradingCurve() const;
	double GetRelativeGrowthRate() const;
	double GetInitialTimeBarrier() const;
	long GetNumEventsMax() const;
	double GetRandomVelocityRange() const;
	double GetTimeMax() const;
	double GetTimePrintOut() const;
	double GetVolumeFraction() const;
	int GetNumThreads() const;
	const std::string& GetDirectory() const;
	void Close();

	bool Is2D() const;
	double GetShrinkage() const;
	double GetAbsoluteGrowthRate() const;


	void OpenFile(std::string rFileName);
	void SkipToNextData();
	double ReadNumber();
	NuTo::FullVector<double, Eigen::Dynamic> ReadVector();
	bool ReadBool();
	std::string ReadString();

	void ReadFile();
	void ReadSimulationParameters();
	void ReadBoundingBox();
	void ReadGradingCurve();

private:

	void CheckInputs() const;
	void mThrow(const std::string& rMsg) const;

	std::ifstream mFile;

	std::string mDirectory;

	long mNumEventsMax;
	double mTimeMax;
	double mTimePrintOut;
	double mInitialTimeBarrier;
	double mRandomVelocityRange;
	double mRelativeGrowthRate;
	double mAbsoluteGrowthRate;
	int mNumThreads;

	int mTypeOfSpecimen;
	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mBoundingBox;
	bool mIs2D;

	NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mGradingCurve;
	double mVolumeFraction;
	double mAbsoluteDistance;

	double mShrinkage;

};

} /* namespace NuTo */
