#pragma once

#include <boost/foreach.hpp>
#include "MoistureTransport_Setup.h"


//! @brief Determines the name of the test
template <int TDim>
std::string SimulationTest_GetTestName(std::array<bool, TDim>& rBS)
{
    std::string testName = std::string("SimulationTest") + std::to_string(TDim) + "D";
    if (TDim > 1)
    {
        testName += "_";
        if (rBS[0])
        {
            testName += "X";
        }
        if (rBS[1])
        {
            testName += "Y";
        }
    }
    if (TDim > 2)
    {
        if (rBS[2])
        {
            testName += "Z";
        }
    }
    return testName;
}


template <int TDim>
void CompareResultsToPaper(NuTo::Structure& rS, std::vector<int> rN, std::vector<double> rL)
{
    int NumActiveDimensions = 0;
    int relevantDirection = 0;
    for (int i = 0; i < TDim; ++i)
    {
        assert((rN[i] == 1 || rN[i] == 16) && "Only 16 or 1 element(s) in each direction allowed for this test");
        if (rN[i] > 1)
        {
            ++NumActiveDimensions;
            relevantDirection = i;
        }
    }
    assert(NumActiveDimensions == 1 && "Results are only valid for one dimensional flows");

    // values fitted from Johannesson and Nyman(2010)
    Eigen::VectorXd PaperValues(17);
    PaperValues[0] = 0.06;
    PaperValues[1] = 0.097;
    PaperValues[2] = 0.116;
    PaperValues[3] = 0.129;
    PaperValues[4] = 0.138;
    PaperValues[5] = 0.146;
    PaperValues[6] = 0.148;
    PaperValues[7] = 0.151;
    PaperValues[8] = 0.152;
    PaperValues[9] = PaperValues[7];
    PaperValues[10] = PaperValues[6];
    PaperValues[11] = PaperValues[5];
    PaperValues[12] = PaperValues[4];
    PaperValues[13] = PaperValues[3];
    PaperValues[14] = PaperValues[2];
    PaperValues[15] = PaperValues[1];
    PaperValues[16] = PaperValues[0];

    typedef boost::ptr_map<int, NuTo::NodeBase> NodeMap;

    assert(rL[relevantDirection] == 0.16 &&
           "The length in flow direction must be 0.16m for direct comparison with paper values");

    double tolerance = 0.005; // Tolerance because not all necessary value (sorption curve) are given in the paper and
    // must be approximated
    double deltaL = rL[relevantDirection] / rN[relevantDirection];
    unsigned int numMismatchingValues = 0;


    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    BOOST_FOREACH (NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;


        if (nodePtr->GetNum(NuTo::Node::eDof::WATERVOLUMEFRACTION) < 1)
        {
            continue; // Nodes without WVF cant be checked --- for example boundary control node
        }

        double relevantNodeCoord = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[relevantDirection];
        int relevantIndex = static_cast<int>(std::round(relevantNodeCoord / deltaL));

        double nodalWVF = nodePtr->Get(NuTo::Node::eDof::WATERVOLUMEFRACTION)[0];
        double paperWVF = PaperValues[relevantIndex];
        if (nodalWVF < paperWVF - tolerance || nodalWVF > paperWVF + tolerance)
        {
            ++numMismatchingValues;
        }
    }
    if (numMismatchingValues > 0)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                              "One ore more calculated values exceeds the tolerance when compared to reference values");
    }
}


//! @brief performs a simulation in the desired dimension
//! @param rN: array with number of elements in each direction
//! @param rL: array with length of elements in each direction
//! @param rBS: array of bools for each direction to determine if the corresponding surfaces are boundary surfaces
//! @param rTestInterpolationTypeCombi if true the simulation does just a single timestep to check if the program works
//! with the current interpolation type combination. Results aren't compared
template <int TDim>
void SimulationTest(std::vector<int> rN, std::vector<double> rL, std::array<bool, TDim> rBS,
                    std::map<NuTo::Node::eDof, NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                    bool rTestInterpolationTypeCombi = false)
{
    std::string testName = SimulationTest_GetTestName<TDim>(rBS);

    if (rTestInterpolationTypeCombi)
    {
        testName += "_TestInterpolationtypeCombi";
    }

    std::cout << std::endl
              << "--------------------------------------------------------------------------" << std::endl
              << "Start test: " << testName << std::endl
              << "--------------------------------------------------------------------------" << std::endl;

    NuTo::Structure S(TDim);
    TimeControl timeControl;

    if (rTestInterpolationTypeCombi)
    {
        testName += "_TestInterpolationtypeCombi";
        timeControl.t_final = 1.0 * 1.0 * 1.0 * 60.0;
        timeControl.t_write = timeControl.t_final * 2;
        timeControl.delta_t = timeControl.t_final;
    }
    else
    {
        timeControl.t_final = 293.0 * 24.0 * 60.0 * 60.0;
        timeControl.t_write = timeControl.t_final;
        timeControl.delta_t = timeControl.t_final / 5.0;
    }


    NuTo::NewmarkDirect TI(&S);


    MoistureTransportControl MT(S);
    MT.InitialRelativeHumidity = 0.95;
    MT.MassExchangeRate = 3.42e-7;
    MT.DiffusionCoefficientRH = 3.9e-10;
    MT.BoundaryEnvironmentalRH = 0.40;


    MT.SetParametersConstitutiveLaw();

    std::string resultDir = std::string("./ConstitutiveLawMoistureTransport_") + testName;

    SetupStructure<TDim>(S, testName);
    auto SEC = SetupSection<TDim>(S);

    auto meshInfo = NuTo::MeshGenerator::Grid(S, rL, rN);

    for (auto it : rDofIPTMap)
        S.InterpolationTypeAdd(meshInfo.second, it.first, it.second);

    S.ElementGroupSetSection(meshInfo.first, SEC);
    S.ElementGroupSetConstitutiveLaw(meshInfo.first, MT.ConstitutiveLawID);

    SetupIntegrationType<TDim>(S, meshInfo.first);

    S.ElementTotalConvertToInterpolationType(); // old used values 1.0e-12,0.001
    MT.ApplyInitialNodalValues();


    auto LambdaGetBoundaryNodes = [rL, rBS](NuTo::NodeBase* rNodePtr) -> bool {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (rBS[0] && ((x >= 0.0 - Tol && x <= 0.0 + Tol) || (x >= rL[0] - Tol && x <= rL[0] + Tol)))
            {
                return true;
            }

            if (TDim > 1)
            {
                double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
                if (rBS[1] && ((y >= 0.0 - Tol && y <= 0.0 + Tol) || (y >= rL[1] - Tol && y <= rL[1] + Tol)))
                {
                    return true;
                }
            }
            if (TDim > 2)
            {
                double z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];
                if (rBS[2] && ((z >= 0.0 - Tol && z <= 0.0 + Tol) || (z >= rL[2] - Tol && z <= rL[2] + Tol)))
                {
                    return true;
                }
            }
        }
        return false;
    }; // GetBoundaryNodesLambda


    auto LambdaTimeDepBoundaryRH = [&MT, &timeControl](double rTime) -> double {
        if (rTime == 0.0)
        {
            return MT.InitialRelativeHumidity;
        }
        else
        {
            if (rTime < timeControl.BC_TransitionTime)
            {
                return MT.InitialRelativeHumidity -
                       sin(rTime / timeControl.BC_TransitionTime * 3.14 / 2.0) *
                               (MT.InitialRelativeHumidity - MT.BoundaryEnvironmentalRH);
            }
            {
                return MT.BoundaryEnvironmentalRH;
            }
        }
    }; // TimeDepBoundaryRHLambda


    SetupConstrainedNodeBoundaryElements<TDim>(S, LambdaGetBoundaryNodes, LambdaTimeDepBoundaryRH);

    MT.SetupStaticData();
    S.NodeBuildGlobalDofs();

    SetupMultiProcessor(S);

    SetupVisualize(S);

    SetupTimeIntegration(TI, timeControl, resultDir);
    NuTo::Timer timer(testName + " --- Timeintegration.Solve");
    TI.Solve(timeControl.t_final);
    if (!rTestInterpolationTypeCombi)
    {
        CompareResultsToPaper<TDim>(S, rN, rL);
    }
}
