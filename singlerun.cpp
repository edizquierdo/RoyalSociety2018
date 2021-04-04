#define PRINTTOFILE
//#define SEED
#define OUTPUT
#define SPEEDOUTPUT
#define COLLIDE

#define ENABLE_CTOR_JSON 1

#include "packages/cxxopts.hpp"

#include "main.h"

int main (int argc, const char* argv[])
{
    double angle = 3.14159 / 2.0;

    if (argc > 1)
    {
        angle = stod(argv[1]);
    }


    cxxopts::Options options("PhysWormSim", "Mechanical and electrophysiological simulation of C. elegans nematode");

    options.add_options()
        ("p,params", "params json file", cxxopts::value<std::string>()->default_value("input/params.json"))
        ("c,coll", "collision tsv file", cxxopts::value<std::string>()->default_value("input/collision_objs.tsv"))
        ("a,angle", "starting angle", cxxopts::value<double>()->default_value(3.14159 / 2.0))
        ("s,seed", "random initialization seen", cxxopts::value<int>()->default_value(0))
        ("o,output", "output dir", cxxopts::value<string>()->default_value("data/run/"))
        ("h,help", "Print usage")
    ;

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
  
    RandomState rs;
    // long seed = static_cast<long>(time(NULL));
    long seed = 0;
    rs.SetRandomSeed(seed);

    std::cout << std::setprecision(10);

    // Code to run simulation:
    InitializeBodyConstants();

    ifstream BestIndividualFile;
    TVector<double> bestVector(1, VectSize);
    BestIndividualFile.open("data/best.gen.dat");
    BestIndividualFile >> bestVector;

    EvaluationFunctionB(bestVector, rs, angle);
    return 0;
}

