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
    BOOST_CHECK_EXCEPTION(timeControl.SetTimestepFunction([]()->double{return 0.0;}),
                          MechanicsException,CheckException);



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

        expectedMsg = "Timestep must be a positive number!";
        BOOST_CHECK_EXCEPTION(timeControl.SetEquidistantTimestepping(0),MechanicsException,CheckException);
        BOOST_CHECK_EXCEPTION(timeControl.SetEquidistantTimestepping(-1337.0),MechanicsException,CheckException);

        timeControl.SetEquidistantTimestepping(1);
        timeControl.Proceed();
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.Proceed(), 3.0);

        timeControl.SetEquidistantTimestepping(5.5);

        timeControl.Proceed();
        timeControl.Proceed();
        BOOST_CHECK_EQUAL(timeControl.Proceed(), 19.5);
    }
    {
        TimeControl timeControl;

//        std::function<double()> myTimeStepFcn = [&timeControl]()->double
//                                        {
//                                            return timeControl.GetCurrentTime();
//                                        };
//        timeControl.SetTimestepFunction(myTimeStepFcn);
    }

}
