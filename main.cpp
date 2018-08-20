// =============================================================
// Evolution of Integrated Neuromechanical Forward Locomotion
// Eduardo Izquierdo
// Indiana University
// February, 2018
// =============================================================

#include <iostream>
#include <iomanip>  // cout precision
#include <math.h>
#include "TSearch.h"
#include "VectorMatrix.h"
#include "Worm.h"

#define EVOLVE
#define PRINTTOFILE
//#define SEED
//#define OUTPUT
//#define SPEEDOUTPUT

using namespace std;

// Integration parameters
const double Duration = 50.0;           // Seconds
const double Transient = 10.0;          //
const double StepSize = 0.01;
const int N_curvs = 23;                 // Number of cuvature points

// Used for Dumping: Frame rate for recording datais set to 50 frames per second
const double fps = 25.0;
const int skip = (int) (1/(StepSize*fps));

// Genotype -> Phenotype Mapping (Ventral cord)
const double	  BiasRange			         	= 15.0;
const double    SCRange                 = 15.0;
const double    CSRange                 = 15.0;
const double    TauMin                  = 0.5; //
const double    TauMax                  = 2.0;

const double    ESRange                 = 2.0;

const double    SRmax                   = 200.0;
const double    NMJmax                  = 1.0;

// (Head)
const double    HCSRange                = 15.0;

// Fitness
const double    AvgSpeed = 0.00022;              // Average speed of the worm in meters per seconds
const double    BBCfit = AvgSpeed*Duration;

// Size of genotype (VC)
int	VectSize = 30;

// ------------------------------------
// Genotype-Phenotype Mapping
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
    // --------------------------------
    // Parameters for the Ventral Nerve Cord Unit
    // --------------------------------
    // Bias
    phen(1) = MapSearchParameter(gen(1), -BiasRange, BiasRange);        // DB, VBa, VBp
    phen(2) = MapSearchParameter(gen(2), -BiasRange, BiasRange);        // DD, VDa, VDp

    // Time Constant
    phen(3) = MapSearchParameter(gen(3), TauMin, TauMax);               // DB, VBa, VBp
    phen(4) = MapSearchParameter(gen(4), TauMin, TauMax);               // DD, VDa, VDp

    // Self connections
    phen(5) = MapSearchParameter(gen(5), -SCRange, SCRange);            // DB, VBa, VBp
    phen(6) = MapSearchParameter(gen(6), -SCRange, SCRange);            // DD, VDa, VDp

    // Chemical synapses
    phen(7) = MapSearchParameter(gen(7), -CSRange, CSRange);            // DB -> DD, VBa -> VDa, VBp -> VDp

    phen(8) = MapSearchParameter(gen(8), -CSRange, CSRange);          // DB -> VDa, DB -> VDp, VBa -> DD /2, VBp -> DD /2

    phen(9) = MapSearchParameter(gen(9), -CSRange, CSRange);          // DD -> VDa

    // Gap junctions across class within unit
    phen(10) = MapSearchParameter(gen(10), 0.0, ESRange);      // DD - VDa, DD - VDp

    // Gap junctions per class
    phen(11) = MapSearchParameter(gen(11), 0.0, ESRange);      // VD - VD, DD - DD
    phen(12) = MapSearchParameter(gen(12), 0.0, ESRange);      // VB - VB, DB - DB

    // Gap junctions across class, across neural unit
    phen(13) = MapSearchParameter(gen(13), 0.0, ESRange);      // VB -> DB+1

    // Stretch receptor
    phen(14) = MapSearchParameter(gen(14), -SRmax, 0.0);        // B- class SR weight

    // NMJ Weight
    phen(15) = MapSearchParameter(gen(15), 0.0, NMJmax);       // DB, VBa, VBp
    phen(16) = MapSearchParameter(gen(16), -NMJmax, 0.0);      // DD, VDa, VDp

    // --------------------------------
    // Parameters for the Head circuit
    // --------------------------------
    // Bias
    phen(17) = MapSearchParameter(gen(17), -BiasRange, BiasRange);    // SMDD, SMDV
    phen(18) = MapSearchParameter(gen(18), -BiasRange, BiasRange);    // RMDD, RMDV

    // Time Constant
    phen(19) = MapSearchParameter(gen(19), TauMin, TauMax);           // SMDD, SMDV
    phen(20) = MapSearchParameter(gen(20), TauMin, TauMax);           // RMDD, RMDV

    // Self connections
    phen(21) = MapSearchParameter(gen(21), -SCRange, SCRange);      // SMDD, SMDV
    phen(22) = MapSearchParameter(gen(22), 4.0, SCRange);           // RMDD, RMDV

    // Chemical synapses
    phen(23) = MapSearchParameter(gen(23), -HCSRange, HCSRange);      // SMDD -> SMDV, SMDV -> SMDD
    phen(24) = MapSearchParameter(gen(24), -HCSRange, HCSRange);      // SMDD -> RMDV, SMDV -> RMDD
    phen(25) = MapSearchParameter(gen(25), -HCSRange, HCSRange);      // RMDD -> RMDV, RMDV -> RMDD

    // Gap junctions across class within unit
    phen(26) = MapSearchParameter(gen(26), 0.0, ESRange);      // SMDD - RMDD, SMDV - RMDV
    phen(27) = MapSearchParameter(gen(27), 0.0, ESRange);      // RMDV - RMDD

    // SMD Stretch Receptor
    phen(28) = MapSearchParameter(gen(28), -SRmax, 0.0);        // SMD- class SR weight

    // NMJ Weight
    phen(29) = MapSearchParameter(gen(29), 0.0, NMJmax);    // SMDD, SMDV
    phen(30) = MapSearchParameter(gen(30), 0.0, NMJmax);    // RMDD, RMDV
}

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

// ------------------------------------
// Display functions
// ------------------------------------
void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
    cout << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}

void ResultsDisplay(TSearch &s)
{
    TVector<double> bestVector;
    ofstream BestIndividualFile;

    bestVector = s.BestIndividual();
    BestIndividualFile.open("best.gen.dat");
    BestIndividualFile << setprecision(32);
    BestIndividualFile << bestVector << endl;
    BestIndividualFile.close();
}

// ------------------------------------
// The main program
// ------------------------------------
#ifdef EVOLVE
int main (int argc, const char* argv[])
{
    std::cout << std::setprecision(10);

    long randomseed = static_cast<long>(time(NULL));

    if (argc == 2)
        randomseed += atoi(argv[1]);

    TSearch s(VectSize);
    TVector<double> phenotype(1, VectSize);

    // save the seed to a file
    ofstream seedfile;
    seedfile.open ("seed.dat");
    seedfile << randomseed << endl;
    seedfile.close();

    // configure the search
    s.SetRandomSeed(randomseed);
    s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
    s.SetSearchResultsDisplayFunction(ResultsDisplay);
    s.SetSelectionMode(RANK_BASED);             //{FITNESS_PROPORTIONATE,RANK_BASED}
    s.SetReproductionMode(GENETIC_ALGORITHM);	// {HILL_CLIMBING, GENETIC_ALGORITHM}
    s.SetPopulationSize(96);
    s.SetMaxGenerations(1000);
    s.SetMutationVariance(0.1);                // For 71 parameters, an estimated avg change of 0.25 for weights (mapped to 15).
    s.SetCrossoverProbability(0.5);
    s.SetCrossoverMode(UNIFORM);              //{UNIFORM, TWO_POINT}
    s.SetMaxExpectedOffspring(1.1);
    s.SetElitistFraction(0.04);
    s.SetSearchConstraint(1);
    s.SetCheckpointInterval(0);
    s.SetReEvaluationFlag(1);

    // redirect standard output to a file
#ifdef PRINTTOFILE
    ofstream evolfile;
    evolfile.open ("fitness.dat");
    cout.rdbuf(evolfile.rdbuf());
#endif

    // Code to run simulation:
    InitializeBodyConstants();

#ifdef SEED
    ifstream BestIndividualFile;
    TVector<double> bestVector(1, VectSize);
    BestIndividualFile.open("best.gen.dat");
    BestIndividualFile >> bestVector;
    s.InitializeSearch();
    for (int i = 1; i <= s.PopulationSize(); i++){
        for (int j = 1; j <= VectSize; j++)
        {
            s.Individual(i)[j] = bestVector[j];
        }
    }
#endif

    s.SetSearchTerminationFunction(NULL);
    s.SetEvaluationFunction(EvaluationFunctionB);
    s.ExecuteSearch();

#ifdef PRINTTOFILE
    evolfile.close();
#endif

    return 0;
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
