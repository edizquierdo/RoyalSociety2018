// =============================================================
// Evolution of Integrated Neuromechanical Forward Locomotion
// Eduardo Izquierdo
// Indiana University
// February, 2018
// 
// Updated by Michael Ivanitskiy (mivanits@umich.edu)
// February, 2021
// =============================================================

#include <iostream>
#include <iomanip>  // cout precision
#include <math.h>
#include "modules/TSearch.h"
#include "modules/VectorMatrix.h"
#include "modules/Worm.h"

using namespace std;

// Integration parameters
const double Duration = 50.0;           // Seconds
const double Transient = 10.0;          //
const double StepSize = 0.01;
const int N_curvs = 23;                 // Number of cuvature points

// Used for Dumping: Frame rate for recording datais set to 50 frames per second
const double fps = 25.0;
const int skip = (int) (1/(StepSize*fps));

void curvRatio(TVector<double> &v, TVector<double> &antposcurv)
{
    for (int i = 1; i <= N_curvs; i++)
    {
        if (i <= 11)
            antposcurv(1) += fabs(v(i));
        else
            antposcurv(2) += fabs(v(i));
    }
}

double EvaluationFunctionB(TVector<double> &v, RandomState &rs)
{
    double fitness;

#ifdef SPEEDOUTPUT
    ofstream fitfile;
    fitfile.open("speed.dat");
#endif

#ifdef OUTPUT
    ofstream bodyfile, actfile, curvfile, paramsfile, voltagefile;
    bodyfile.open("body.dat");
    actfile.open("act.dat");
    curvfile.open("curv.dat");
    paramsfile.open("params.dat");
#endif

    // Fitness
    fitness = 0.0;
    double bodyorientation, anglediff;
    double movementorientation, distancetravelled = 0, temp;
    TVector<double> curvature(1, N_curvs);
    TVector<double> antpostcurv(1, 2);
    antpostcurv.FillContents(0.0);

    // Genotype-Phenotype Mapping
    TVector<double> phenotype(1, VectSize);
    GenPhenMapping(v, phenotype);

    Worm w(phenotype, 0);

#ifdef OUTPUT
    w.DumpParams(paramsfile);
#endif

    w.InitializeState(rs);

    // Transient
    for (double t = 0.0; t <= Transient; t += StepSize)
    {
        w.Step(StepSize, 1);
#ifdef OUTPUT
        w.Curvature(curvature);
        curvfile << curvature << endl;
        w.DumpBodyState(bodyfile, skip);
        w.DumpActState(actfile, skip);
#endif
    }

    double xt = w.CoMx(), xtp;
    double yt = w.CoMy(), ytp;

    // Time loop
    for (double t = 0.0; t <= Duration; t += StepSize) {

        w.Step(StepSize, 1);

        // Current and past centroid position
        xtp = xt; ytp = yt;
        xt = w.CoMx(); yt = w.CoMy();

        // Integration error check
        if (isnan(xt) || isnan(yt) || sqrt(pow(xt-xtp,2)+pow(yt-ytp,2)) > 100*AvgSpeed*StepSize)
        {
            return 0.0;
        }

        // Fitness
        bodyorientation = w.Orientation();                  // Orientation of the body position
        movementorientation = atan2(yt-ytp,xt-xtp);         // Orientation of the movement
        anglediff = movementorientation - bodyorientation;  // Check how orientations align
        temp = cos(anglediff) > 0.0 ? 1.0 : -1.0;           // Add to fitness only movement forward
        distancetravelled += temp * sqrt(pow(xt-xtp,2)+pow(yt-ytp,2));

#ifdef OUTPUT
        w.Curvature(curvature);
        curvfile << curvature << endl;
        w.DumpBodyState(bodyfile, skip);
        w.DumpActState(actfile, skip);
#endif
    }
    fitness = 1 - (fabs(BBCfit-distancetravelled)/BBCfit);

#ifdef OUTPUT
    cout << fitness << " " << BBCfit << " " << distancetravelled << " " << distancetravelled/Duration << endl;
    bodyfile.close();
    actfile.close();
    curvfile.close();
#endif

#ifdef SPEEDOUTPUT
    fitfile << fitness << " "<< BBCfit << " " << distancetravelled << " " << distancetravelled/Duration << " " << endl;
    fitfile.close();
#endif

    return fitness;
}

#else
int main (int argc, const char* argv[])
{
    RandomState rs;
    long seed = static_cast<long>(time(NULL));
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
#endif
