#include "BoostUnitTest.h"

#include "mechanics/timeIntegration/TimeControl.h"

#include <iostream>

using namespace NuTo;

// %%% Helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool CheckExceptionMessageCorrect(MechanicsException const& ex,std::string expectedMsg)
{
    std::string refinedErrorMessage = ex.ErrorMessage().substr(ex.ErrorMessage().find_first_of("]")+2);

    if (expectedMsg.compare(refinedErrorMessage))
    {
        std::cout << std::endl
                  << "WRONG EXCEPTION TRIGGERED!" << std::endl
                  << "Expected exception message:" << std::endl << expectedMsg << std::endl
                  << "Triggered exception message:" << std::endl << refinedErrorMessage << std::endl << std::endl;
        return false;
    }
    return true;
}


// %%% Tests %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BOOST_AUTO_TEST_CASE(Debug_Fuctionality)
{
    TimeControl timeControl;

    std::string expectedMsg = "";
    std::function<bool(MechanicsException const&)> CheckException
            = [&expectedMsg](MechanicsException const& ex)->bool
                {
                    return CheckExceptionMessageCorrect(ex,expectedMsg);
                };

    expectedMsg = "Current timestep is 0 or negative!";
//    BOOST_CHECK_EXCEPTION(timeControl.SetTimeStepFunction([](double,double)->double{return 0.0;}),
//                          MechanicsException,CheckException);



}






BOOST_AUTO_TEST_CASE(Equidistant_Timestepping)
{
    std::string expectedMsg = "";
    std::function<bool(MechanicsException const&)> CheckException
            = [&expectedMsg](MechanicsException const& ex)->bool
                {
                    return CheckExceptionMessageCorrect(ex,expectedMsg);
                };
    {
        TimeControl timeControl;
        timeControl.SetTimeFinal(20);

        expectedMsg = "Timestep must be larger tham 0!";
        BOOST_CHECK_EXCEPTION(timeControl.SetTimeStep(0),MechanicsException,CheckException);
        BOOST_CHECK_EXCEPTION(timeControl.SetTimeStep(-1337.0),MechanicsException,CheckException);

        timeControl.SetTimeStep(1);
        timeControl.Proceed();
        timeControl.Proceed();
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.GetCurrentTime(), 3.0);

        timeControl.SetTimeStep(5.5);

        timeControl.Proceed();
        timeControl.Proceed();
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.GetCurrentTime(), 19.5);
    }
    {
        TimeControl timeControl;
        timeControl.SetTimeFinal(20);
//        std::function<double()> myTimeStepFcn = [&timeControl]()->double
//                                        {
//                                            return timeControl.GetCurrentTime();
//                                        };
//        timeControl.SetTimestepFunction(myTimeStepFcn);
    }

}
