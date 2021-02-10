// *****************************
// "Evolutionary" search classes
// *****************************

// Uncomment the following line to enable multithreading
#define THREADED_SEARCH
#define THREAD_COUNT 8

#pragma once

#include "VectorMatrix.h"
#include "random.h"
#ifdef THREADED_SEARCH
  #include <pthread.h>
#endif

using namespace std;


// A utility function for clipping a double to lie within the interval [MIN,MAX]

inline double clip(double x, double min, double max)
{
	double temp;

	temp = ((x > min)?x:min);
	return (temp < max)?temp:max;
}


// The minimum and maximum allowable parameter values

const double MinSearchValue = -1.0, MaxSearchValue = 1.0;


// A utility function for mapping from search parameters to model parameters,
// with the result optionally clipped to the given range.

inline double MapSearchParameter(double x, double min, double max,
                                           double clipmin = -1.0e99, double clipmax = 1.0e99)
{
	double m = (max - min)/(MaxSearchValue - MinSearchValue);
	double b = min - m * MinSearchValue;

	return clip(m * x + b,clipmin,clipmax);
}


// A utility function for mapping from model parameters to search parameters

inline double InverseMapSearchParameter(double x, double min, double max)
{
	double m = (MaxSearchValue - MinSearchValue)/(max - min);
	double b = MinSearchValue - m * min;

	return m * x + b;
}


// *******************************
// The TSearch class declaration
// *******************************

enum TSelectionMode {FITNESS_PROPORTIONATE,RANK_BASED};    // Supported selection modes
enum TReproductionMode {HILL_CLIMBING, GENETIC_ALGORITHM}; // Supported reproduction modes
enum TCrossoverMode {UNIFORM, TWO_POINT};                  // Supported crossover modes

class TSearch {
	public:
		// The constructor
		TSearch(int vectorSize = 0,double (*EvalFn)(TVector<double> &, RandomState &) = NULL);
		// The destructor
		~TSearch();
		// Basic Accessors
		int VectorSize(void) {return vectorSize;};
		void SetVectorSize(int NewSize);
        void SetRandomSeed(long seed) {rs.SetRandomSeed(seed);};
		// Search Mode Accessors
		TSelectionMode SelectionMode(void) {return SelectMode;};
		void SetSelectionMode(TSelectionMode newmode) {SelectMode = newmode;};
		TReproductionMode ReproductionMode(void) {return RepMode;};
		void SetReproductionMode(TReproductionMode newmode) {RepMode = newmode;};
		TCrossoverMode CrossoverMode(void) {return CrossMode;};
		void SetCrossoverMode(TCrossoverMode newmode) {CrossMode = newmode;};
		// Search Parameter Accessors
		int PopulationSize(void) {return Population.Size();};
		void SetPopulationSize(int NewSize);
		int MaxGenerations(void) {return MaxGens;};
		void SetMaxGenerations(int NewMax);
		double ElitistFraction(void) {return EFraction;};
		void SetElitistFraction(double NewFraction);
		double MaxExpectedOffspring(void) {return MaxExpOffspring;};
		void SetMaxExpectedOffspring(double NewVal);
		double MutationVariance(void) {return MutationVar;};
		void SetMutationVariance(double NewVariance);
		double CrossoverProbability(void) {return CrossProb;};
		void SetCrossoverProbability(double NewProb);
		TVector<int> &CrossoverTemplate(void) {return crossTemplate;};
		void SetCrossoverTemplate(TVector<int> &NewTemplate);
		TVector<int> &CrossoverPoints(void) {return crossPoints;};
		void SetCrossoverPoints(TVector<int> &NewPoints);
		TVector<int> &SearchConstraint(void) {return ConstraintVector;};
		void SetSearchConstraint(TVector<int> &Constraint);
		void SetSearchConstraint(int Flag);
		int ReEvaluationFlag(void) {return ReEvalFlag;};
		void SetReEvaluationFlag(int flag) {ReEvalFlag = flag;};
		double CheckpointInterval(void) {return CheckpointInt;};
		void SetCheckpointInterval(int NewFreq);
		// Function Pointer Accessors
		void SetEvaluationFunction(double (*EvalFn)(TVector<double> &v, RandomState &rs))
			{EvaluationFunction = EvalFn;};
		void SetBestActionFunction(void (*BestFn)(int Generation,TVector<double> &v))
			{BestActionFunction = BestFn;};
		void SetPopulationStatisticsDisplayFunction(void (*DisplayFn)(int Generation,double BestPerf,double AvgPerf,double PerfVar))
			{PopulationStatisticsDisplayFunction = DisplayFn;};
		void SetSearchTerminationFunction(int (*TerminationFn)(int Generation,double BestPerf,double AvgPerf,double PerfVar))
			{SearchTerminationFunction = TerminationFn;};
		void SetSearchResultsDisplayFunction(void (*DisplayFn)(TSearch &s))
			{SearchResultsDisplayFunction = DisplayFn;};
		// Status Accessors
		int Generation(void) {return Gen;};
		TVector<double> &Individual(int i) {return Population(i);};
		double Fitness(int i) {return fitness(i);};
		double Performance(int i) {return Perf(i);};
		double BestPerformance (void) {return BestPerf;};
		TVector<double> &BestIndividual(void) {return bestVector;};
		// Control
		void InitializeSearch(void);
		void ExecuteSearch(void);
		void ResumeSearch(void);
		// Input and output
    void WriteCheckpointFile(void);
    void ReadCheckpointFile(void);
    //friend ostream& operator<<(ostream& os, TSearch& s);
		//friend istream& operator>>(istream& is, TSearch& s);

	private:
		// Helper Methods
		void DoSearch(int ResumeFlag);
		int EqualVector(TVector<double> &v1, TVector<double> &v2)
		{
			if (v1.Size() != v2.Size()) return 0;
			if (v1.LowerBound() != v2.LowerBound()) return 0;
			if (v1.UpperBound() != v2.UpperBound()) return 0;
			for (int i = v1.LowerBound(); i <= v1.UpperBound(); i++)
				if (v1[i] != v2[i]) return 0;
			return 1;
		};
		void RandomizeVector(TVector<double> &Vector);
		void RandomizePopulation(void);
		double EvaluateVector(TVector<double> &Vector, RandomState &rs);
    friend void *EvaluatePopulationRange(void *arg);
		void EvaluatePopulation(int start = 1);
		void SortPopulation(void);
		void UpdatePopulationFitness(void);
		void ReproducePopulationHillClimbing(void);
		void ReproducePopulationGeneticAlgorithm(void);
		void MutateVector(TVector<double> &Vector);
		void UniformCrossover(TVector<double> &v1, TVector<double> &v2);
		void TwoPointCrossover(TVector<double> &v1, TVector<double> &v2);
		void PrintPopulationStatistics(void);
		void ReproducePopulation(void);
		void UpdatePopulationStatistics(void);
		void DisplayPopulationStatistics(void);
		int SearchTerminated(void);
		void DisplaySearchResults(void);

		// Internal State
    RandomState rs;
    TVector<RandomState> RandomStates;
		int Gen;
		int SearchInitialized;
		TVector<TVector<double> > Population;
		TVector<double> Perf;
		TVector<double> fitness;
		int UpdateBestFlag;
		TVector<double> bestVector;
		double BestPerf;
		double MinPerf, MaxPerf, AvgPerf, PerfVar;
		// Search Modes
		TSelectionMode SelectMode;
		TReproductionMode RepMode;
		TCrossoverMode CrossMode;
		// Search Parameters
		int vectorSize;
		int MaxGens;
		double EFraction;
		double MaxExpOffspring;
		double MutationVar;
		double CrossProb;
		TVector<int> crossTemplate;
		TVector<int> crossPoints;
		TVector<int> ConstraintVector;
		int ReEvalFlag;
		int CheckpointInt;
		// Function Pointers
		double (*EvaluationFunction)(TVector<double> &v, RandomState &rs);
		void (*BestActionFunction)(int Generation,TVector<double> &v);
		void (*PopulationStatisticsDisplayFunction)(int Generation,double BestPerf,double AvgPerf,double PerfVar);
		int (*SearchTerminationFunction)(int Generation,double BestPerf,double AvgPerf,double PerfVar);
		void (*SearchResultsDisplayFunction)(TSearch &s);
};


// A range specification structure for threaded evaluation

struct PopRangeSpec {TSearch *search; int start; int end;};
