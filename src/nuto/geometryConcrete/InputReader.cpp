/*
 * InputReader.cpp
 *
 *  Created on: 17 Mar 2014
 *      Author: ttitsche
 */

#include "nuto/geometryConcrete/InputReader.h"
#include "nuto/base/Exception.h"
#include <sstream>

NuTo::InputReader::InputReader(std::string rFileName)
{
	OpenFile(rFileName);

	ReadFile();

	mFile.close();
}

void NuTo::InputReader::OpenFile(std::string rFileName)
{
	mFile.open(rFileName.c_str(), std::ios::in);
	if (!mFile.is_open())
		throw Exception("[NuTo::InputReader::OpenFile] File " "" + rFileName + "" " not found.");
}

void NuTo::InputReader::ReadFile()
{
	ReadSimulationParameters();
	ReadBoundingBox();

	mIsDefinedByGradingCurve = ReadBool();

	if (mIsDefinedByGradingCurve)
		ReadGradingCurve();
	else
		ReadMonodisperse();

	mShrinkage = ReadNumber();
}

void NuTo::InputReader::ReadSimulationParameters()
{
	mDirectory = ReadString();

	mNumEventsMax = ReadNumber();
	mTimeMax = ReadNumber();
	mTimePrintOut = ReadNumber();
	mInitialTimeBarrier = ReadNumber();
	mRandomVelocityRange = ReadNumber();
	mGrowthRates = ReadNumber();
	mNumThreads = ReadNumber();
}

void NuTo::InputReader::ReadBoundingBox()
{
	mBoxType = ReadNumber();
	FullVector<double, Eigen::Dynamic> boxVector = ReadVector();
	if (boxVector.GetNumRows() != 6)
	{
		std::stringstream exceptionStream;
		exceptionStream << "[NuTo::InputReader::ReadFile] boundingBox - " << boxVector.GetNumRows() << " components in input file, expected 6.";
		throw Exception(exceptionStream.str());
	}
	mBoundingBox = FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>(3, 2);
	mBoundingBox << boxVector[0], boxVector[1], boxVector[2], boxVector[3], boxVector[4], boxVector[5];
}

void NuTo::InputReader::ReadGradingCurve()
{
	int nSieves = ReadNumber();
	mGradingCurve = FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>(nSieves, 3);
	for (int i = 0; i < nSieves; ++i)
	{
		FullVector<double, Eigen::Dynamic> gradingCurveVector = ReadVector();
		if (gradingCurveVector.GetNumRows() != 3)
		{
			std::stringstream exceptionStream;
			exceptionStream << "[NuTo::InputReader::ReadVector] gradingCurve - " << gradingCurveVector.GetNumRows() << " components in input file, expected 3.";
			throw Exception(exceptionStream.str());
		}
		mGradingCurve.SetRow(i, gradingCurveVector.transpose());
	}
	mVolumeFraction = ReadNumber();
	mAbsoluteDistance = ReadNumber();

}

void NuTo::InputReader::ReadMonodisperse()
{
	mNumParticles = ReadNumber();
}

void NuTo::InputReader::SkipToNextData()
{
	const int skipLength = 1000;
	char skipText[skipLength];
	char skipEqualSign;
	mFile.get(skipText, skipLength, '=');
	mFile.get(skipEqualSign);
}

bool NuTo::InputReader::ReadBool()
{
	SkipToNextData();


	char value;
	std::vector<char> positive = { 'T', 't', 'Y', 'y', '1' };
	mFile >> value;
	for (char& c : positive)
		if (c == value)
			return true;

	return false;
}

double NuTo::InputReader::ReadNumber()
{
	SkipToNextData();

	double number;
	mFile >> number;
	return number;
}

NuTo::FullVector<double, Eigen::Dynamic> NuTo::InputReader::ReadVector()
{
	SkipToNextData();

	std::vector<double> readVector;
	double number;

	char skipComma = ',';

	while (skipComma == ',')
	{
		mFile >> number;
		readVector.push_back(number);
		mFile.get(skipComma);
	}

	return readVector;
}

std::string NuTo::InputReader::ReadString()
{
	SkipToNextData();
	char directory[1000];
	mFile.get(directory,';');

	std::string str = directory;

	// remove ';'
	str.erase(str.size()-1);

	// remove whitespaces
	str.erase(remove_if(str.begin(), str.end(), isspace), str.end());

	return str;
}

void NuTo::InputReader::PrintInput()
{

	std::cout << "          INPUT PARAMETERS" << std::endl;
	std::cout << "======================================" << std::endl;
	std::cout << std::endl <<  "     SIMULATION PARAMETERS   " << std::endl;
	std::cout << " directory:                  " << mDirectory << std::endl;
	std::cout << " max. number of events:      " << mNumEventsMax << std::endl;
	std::cout << " max. simulation time [s]:   " << mTimeMax << std::endl;
	std::cout << " time between printouts [s]: " << mTimePrintOut << std::endl;
	std::cout << " initial time barrier:       " << mInitialTimeBarrier << std::endl;
	std::cout << " random velocity range:      " << mRandomVelocityRange << std::endl;
	std::cout << " growth rates:               " << mGrowthRates << std::endl;
	std::cout << std::endl <<  "        BOUNDING BOX" << std::endl;
	std::cout << "    box type:                " << mBoxType << std::endl;
	std::cout << "    x-Range  " << mBoundingBox.GetValue(0, 0) << " to " << mBoundingBox.GetValue(0, 1) << std::endl;
	std::cout << "    y-Range  " << mBoundingBox.GetValue(1, 0) << " to " << mBoundingBox.GetValue(1, 1) << std::endl;
	std::cout << "    z-Range  " << mBoundingBox.GetValue(2, 0) << " to " << mBoundingBox.GetValue(2, 1) << std::endl;
	std::cout << std::endl << "          PARTICLES   " << std::endl;
	if (mIsDefinedByGradingCurve)
	{
		std::cout << "The particles are defined by the following grading curve: " << std::endl;
		std::cout << mGradingCurve << std::endl;
		std::cout << " volume fraction: " << mVolumeFraction << std::endl;
		std::cout << " absolute distance: " << mAbsoluteDistance << std::endl;
	}
	else
	{
		std::cout << "The " << mNumParticles << " particles are initialized with radius = 0." << std::endl;
	}

}

double NuTo::InputReader::GetAbsoluteDistance() const
{
	return mAbsoluteDistance;
}

const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& NuTo::InputReader::GetBoundingBox() const
{
	return mBoundingBox;
}

int NuTo::InputReader::GetBoxType() const
{
	return mBoxType;
}

const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& NuTo::InputReader::GetGradingCurve() const
{
	return mGradingCurve;
}

double NuTo::InputReader::GetGrowthRates() const
{
	return mGrowthRates;
}

double NuTo::InputReader::GetInitialTimeBarrier() const
{
	return mInitialTimeBarrier;
}

bool NuTo::InputReader::isIsDefinedByGradingCurve() const
{
	return mIsDefinedByGradingCurve;
}

long NuTo::InputReader::GetNumEventsMax() const
{
	return mNumEventsMax;
}

int NuTo::InputReader::GetNumParticles() const
{
	return mNumParticles;
}

double NuTo::InputReader::GetRandomVelocityRange() const
{
	return mRandomVelocityRange;
}

double NuTo::InputReader::GetTimeMax() const
{
	return mTimeMax;
}

double NuTo::InputReader::GetTimePrintOut() const
{
	return mTimePrintOut;
}

double NuTo::InputReader::GetVolumeFraction() const
{
	return mVolumeFraction;
}

int NuTo::InputReader::GetNumThreads() const
{
	return mNumThreads;
}

const std::string& NuTo::InputReader::GetDirectory() const
{
	return mDirectory;
}
