#include "LinearElasticBenchmarkStructure.h"
#include "base/Timer.h"

int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        std::cout << "Please provide the following arguments for this benchmark\n"
                  << "1) number of elements in x direction\n"
                  << "2) number of elements in x direction\n"
                  << "3) number of elements in x direction\n"
                  << "4) number of OPENMP processors\n";
        return 0;
    }
    const int numX = std::stoi(argv[1]);
    const int numY = std::stoi(argv[2]);
    const int numZ = std::stoi(argv[3]);
    const int numProc = std::stoi(argv[4]);

    NuTo::Timer timerTotal("TotalTime", true);

    NuTo::Timer timer("LinearElasticityCustom::Setup", true);
    NuTo::Benchmark::LinearElasticBenchmarkStructure s({numX, numY, numZ}, numProc);

    timer.Reset("LinearElasticityCustom::Newmark");
    s.SolveWithNewmark(1, "LinearElasticCustomResults");
    std::cout << " #### Structure::GetNumActiveDofs() " << s.GetStructure().GetNumTotalActiveDofs() << " ####\n";

    return 0;
}
