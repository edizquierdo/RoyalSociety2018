#define PRINTTOFILE
//#define SEED
#define OUTPUT
#define SPEEDOUTPUT
#define COLLIDE

#define ENABLE_CTOR_JSON 1

#include "main.h"

int main (int argc, const char* argv[])
{
    double angle = 3.14159 / 2.0;

    if (argc > 1)
    {
        angle = stod(argv[1]);
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

