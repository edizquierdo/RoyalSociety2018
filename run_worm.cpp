// #define PRINTTOFILE
//#define SEED
//#define OUTPUT
//#define SPEEDOUTPUT

#include "main.h"

int main (int argc, const char* argv[])
{
    RandomState rs;
    // long seed = static_cast<long>(time(NULL));
    long seed = 0;
    rs.SetRandomSeed(seed);

    std::cout << std::setprecision(10);

    // Code to run simulation:
    InitializeBodyConstants();

    ifstream BestIndividualFile;
    TVector<double> bestVector(1, VectSize);
    BestIndividualFile.open("best.gen.dat");
    BestIndividualFile >> bestVector;

    EvaluationFunctionB(bestVector, rs);
    return 0;
}

