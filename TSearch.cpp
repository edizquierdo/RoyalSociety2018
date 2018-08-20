// *******************************************************************************
// Methods for the evolutionary search class TSearch
//
// RDB 
//  1/99 - Created
//  5/07 - Added binary checkpoint files (with contributions from Chad Seys)
//  1/08 - Added multithreaded evaluation (with contributions from Chad Seys and Paul Williams)
//
// TO DO
//   1. Abstract TSearch over the type of individuals, so that more than just
//      real vectors can be searched
//   2. Define specialized support (either a subclass of TSearch or of Individual)
//      for evolving CTRNNs. Features might include support for setting parameter
//      ranges, seeding w/ center-crossing, setting up symmetric circuits, 
//      automating search-to-CTRNN parameter mapping, etc.
//   3. Add support for co-evolution by allowing two search objects to
//      interact with one another during evolution.  Much of this can probably
//      just be handled by making both search objects global and having each
//      evaluation function refer to the population in the other object. 
//      However, the generations of the search objects must also be interleaved.
//      For example, we could have another function that kept reseting MaxGens
//      and calling ExecuteSearch for each object.
// *******************************************************************************

#include "TSearch.h"
#include <math.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>


// An out of memory handler for new

#include <new>

void OutOfMemoryHandler(void)
{
	cerr << "Error: Out of memory!\n";
	exit(0);
}


// *****************************
// Constructors and Destructors
// *****************************

// The constructor

TSearch::TSearch(int VSize, double (*EvalFn)(TVector<double> &, RandomState &))
{
	// Install a default new handler if none is currently installed
	new_handler OldHandler;
	OldHandler = set_new_handler(OutOfMemoryHandler);
	if (OldHandler != NULL) set_new_handler(OldHandler);
	// Initialize internal state                                              
	SearchInitialized = 0;
	// Initialize function pointers
	EvaluationFunction = EvalFn;
	BestActionFunction = NULL;
	SearchTerminationFunction = NULL;
	PopulationStatisticsDisplayFunction = NULL;
	SearchResultsDisplayFunction = NULL;
	// Initialize the vector size
	SetVectorSize(VSize);
	// Set up search mode defaults
	SetSelectionMode(RANK_BASED);
	SetReproductionMode(GENETIC_ALGORITHM);
	SetCrossoverMode(TWO_POINT);
	// Set up search parameter defaults
	SetPopulationSize(1);
	SetMaxGenerations(0);
	SetElitistFraction(0.0);
	SetMaxExpectedOffspring(1.1);
	SetMutationVariance(1.0);
	SetCrossoverProbability(0.0);
	SetSearchConstraint(1);
	SetReEvaluationFlag(0);
	SetCheckpointInterval(0);
}


// The destructor

TSearch::~TSearch()
{
  RandomStates.SetSize(0);
	for (int i = 1; i <= PopulationSize(); i++)
		Population[i].SetSize(0);
	Population.SetSize(0);
	Perf.SetSize(0);
	fitness.SetSize(0);
	crossTemplate.SetSize(0);
	crossPoints.SetSize(0);
	ConstraintVector.SetSize(0);
	bestVector.SetSize(0);
}


// *********
// Accessors
// *********

// Resize the search vector and related vectors

void TSearch::SetVectorSize(int NewSize)
{
	// Set up the new vector size
	if (NewSize <= 0) {cerr << "Invalid vector size: "<< NewSize; exit(0);}
	vectorSize = NewSize;
	// Resize the population
	for (int i = 1; i <= Population.Size(); i++)
		Population[i].SetSize(NewSize);
	// Adjust bestVector
	bestVector.SetSize(NewSize);
	// Reset the crossover template and crossover points vectors
	TVector<int> v(1,NewSize);
	for (int i = 1; i <= NewSize; i++)
		v[i] = i;
	SetCrossoverTemplate(v);
	// Reset the constraint vector
	ConstraintVector.SetSize(NewSize);
	ConstraintVector.FillContents(1);
}


// Resize the population vector and related vectors

void TSearch::SetPopulationSize(int NewSize)
{
	if (NewSize <= 0) {cerr << "Invalid population size: "<< NewSize; exit(0);}
	Population.SetSize(NewSize);
	for (int i = 1; i <= NewSize; i++)
		Population[i].SetSize(vectorSize);
	Perf.SetSize(NewSize);
	fitness.SetSize(NewSize);
  RandomStates.SetSize(NewSize);
  for (int i = 1; i <= NewSize; i++)
      RandomStates[i].SetRandomSeed(rs.UniformRandomInteger(1,32767)); // XXX
//    RandomStates[i].SetRandomSeed(rs.UniformRandomInteger(1,LONG_MAX));
}


// Set the maximum number of generations to search

void TSearch::SetMaxGenerations(int NewMax) 
{
	if (NewMax < 0) {
		cerr << "Invalid MaxGenerations: " << NewMax; 
		exit(0);
	}
	MaxGens = NewMax;
}


// Set the fraction of the new population to be produced via elitist selection

void TSearch::SetElitistFraction(double NewFraction)
{
	if (NewFraction < 0.0 || NewFraction > 1.0) {
		cerr << "Invalid ElitismFraction: " << NewFraction; 
		exit(0);
	}
	EFraction = NewFraction;
}
		

// Set the number of offspring to be allocated to the highest-performing individual

void TSearch::SetMaxExpectedOffspring(double NewVal)
{
	if (NewVal < 1.0 || NewVal > 2.0) {
		cerr << "Invalid MaxExpectedOffspring: " << NewVal; 
		exit(0);
	}
	MaxExpOffspring = NewVal;
}
		

// Set the mutation variance

void TSearch::SetMutationVariance(double NewVariance) 
{
	if (NewVariance <= 0.0) {
		cerr << "Invalid MutationVariance: " << NewVariance; 
		exit(0);
	}
	MutationVar = NewVariance;
}
		

// Set the crossover probability

void TSearch::SetCrossoverProbability(double NewProb) 
{
	if (NewProb < 0.0 || NewProb > 1.0) {
		cerr << "Invalid CrossoverProbability: " << NewProb; 
		exit(0);
	}
	CrossProb = NewProb;
}
		

// Set the crossover template

void TSearch::SetCrossoverTemplate(TVector<int> &NewTemplate)
{
	// Modify CrossoverTemplate
	if (NewTemplate.Size() != vectorSize) {
		cerr << "Invalid vector size for CrossoverTemplate: " << NewTemplate.Size();
		exit(0);
	}
	int x = 1;
	for (int i = 1; i <= NewTemplate.Size(); i++)
        if (NewTemplate[i] != x)
        {
            if (NewTemplate[i] == x+1)
            {
                x++;
            }
			else
            {
				cerr << "Invalid format for CrossoverTemplate: " << NewTemplate;
				exit(0);
			}
        }
	crossTemplate = NewTemplate;
	// Modify CrossoverPoints appropriately
	crossPoints.SetSize(x);
	crossPoints[1] = 1;
	x = 1;
	for (int i = 1; i <= vectorSize; i++)
		if (NewTemplate[i] != x)
			crossPoints[++x] = i;
}


// Set the crossover points

void TSearch::SetCrossoverPoints(TVector<int> &NewPoints)
{
	// Modify CrossoverPoints
	if (NewPoints.Size() > vectorSize) {
		cerr << "Invalid vector size for Crossover Points: " << NewPoints.Size();
		exit(0);
	}
	if (NewPoints.Size() < 1 || NewPoints[1] != 1) {
		cerr << "Invalid format for Crossover Points: " << NewPoints;
		exit(0);
	}
	int x = 0;
	for (int i = 1; i <= NewPoints.Size(); i++)
		if (NewPoints[i] > x && NewPoints[i] <= vectorSize) x = NewPoints[i];
		else {
			cerr << "Invalid format for Crossover Points: " << NewPoints;
			exit(0);
		}
	crossPoints = NewPoints;
	// Modify CrossoverTemplate appropriately
	x = 1;
	for (int i = 1; i < NewPoints.Size(); i++) {
		for (int j = NewPoints[i]; j <= NewPoints[i+1]; j++) {
			int k = (j <= vectorSize)?j:vectorSize;
			crossTemplate[k] = x;
		}
		x++;
	}
	for (int i = NewPoints[NewPoints.Size()]; i <= vectorSize; i++)
		crossTemplate[i] = x;
}


// Set the search constraint vector

void TSearch::SetSearchConstraint(TVector<int> &constraint) 
{
	if (constraint.Size() != vectorSize) {
		cerr << "Invalid vector size for SearchConstraint: " << constraint;
		exit(0);
	}
	ConstraintVector = constraint;
}

void TSearch::SetSearchConstraint(int flag) 
{
	ConstraintVector.FillContents(flag);
}


// Set the frequency with which checkpoint files are written
// (0 means never)

void TSearch::SetCheckpointInterval(int NewInterval) 
{
	if (NewInterval < 0) {
		cerr << "Invalid CheckpointInterval: " << NewInterval; 
		exit(0);
	}
	CheckpointInt = NewInterval;
}


// *****************
// Basic Search Loop
// *****************

// The top-level search loop

void TSearch::DoSearch(int ResumeFlag)
{
	// Initialize search if necessary
	if (!SearchInitialized) InitializeSearch();
	// Make sure we have an evaluation function
	if (EvaluationFunction == NULL)
	{
		cerr << "Error: NULL evaluation function\n";
		exit(0);
	}
	// Unless we're resuming a checkpointed search, evalute the initial population and reset best
	if (!ResumeFlag) {
		EvaluatePopulation();
		BestPerf = -1;
		UpdateBestFlag = 0;
	}
	// Update and display statistics of the initial population
	UpdatePopulationStatistics();
	DisplayPopulationStatistics();
	// If the best changed and there is a BestActionFunction, invoke it
	if (UpdateBestFlag && BestActionFunction != NULL)
		(*BestActionFunction)(Gen,bestVector);
	// Repeat until done
	while (!SearchTerminated())
	{
		Gen++;
		UpdateBestFlag = 0;
		ReproducePopulation();
		UpdatePopulationStatistics();
		DisplayPopulationStatistics();
		// If the best changed and there is a BestActionFunction, invoke it
		if (UpdateBestFlag && BestActionFunction != NULL)
			(*BestActionFunction)(Gen,bestVector);
		// If we're checkpointing and this is a checkpoint generation, save the state of the search 
		if ((CheckpointInt > 0) && (Gen > 0) && ((Gen % CheckpointInt) == 0))
			WriteCheckpointFile();
	}
	// Display results
	DisplaySearchResults();
}


// Execute a search

void TSearch::ExecuteSearch(void)
{
	DoSearch(0);
}


// Initialize a new search

void TSearch::InitializeSearch(void)
{
	// Reset the generation counter
	Gen = 0;
	// Set up the initial population
	RandomizePopulation();
	// The search is now initialized
	SearchInitialized = 1;
}


// Randomize a vector

void TSearch::RandomizeVector(TVector<double> &v)
{
	for (int i = 1; i <= v.Size(); i++)
		v[i] = rs.UniformRandom(MinSearchValue,MaxSearchValue);
}


// Randomize the population

void TSearch::RandomizePopulation(void)
{
	for (int i = 1; i <= Population.Size(); i++)
		RandomizeVector(Population[i]);
}


// Update BestVector, BestPerformance, MinPerformance, MaxPerformance, 
// AveragePerformance and Performance Variance

void TSearch::UpdatePopulationStatistics(void)
{
	register int i;
	double total = 0;
	int bestindex = 1;
	register double perf;

	// Collect various info about the current population
	MinPerf = 1E10;
	MaxPerf = -1E10;
	for (i = 1; i <= Population.Size(); i++)
	{
		perf = Perf[i];
		// Update MinPerformance and MaxPerformance as necessary
		if (perf > MaxPerf) {MaxPerf = perf; bestindex = i;}
		if (perf < MinPerf) MinPerf = perf;
		// Update total
		total += perf;
	}
	// Update AveragePerformance (with protection from possible numerical errors)
	AvgPerf = total/Population.Size();
	if (AvgPerf < MinPerf) AvgPerf = MinPerf;
	if (AvgPerf > MaxPerf) AvgPerf = MaxPerf;
	// Update PerformanceVariance
	if (Population.Size() > 1)
	{ 
		total = 0;
		for (int i = 1; i <= Population.Size(); i++) {
			double d = Perf[i] - AvgPerf;
			total += d*d;
		}
		PerfVar = total/(Population.Size()-1);
	}
	else PerfVar = 0.0;
	// If the best performance has improved or ReEvalFlag is set, update BestPerf and BestVector
	if ((MaxPerf > BestPerf) || ReEvalFlag)
	{
		UpdateBestFlag = 1;
		BestPerf = MaxPerf;
		bestVector = Population[bestindex];
	}
}


// Display population statistics

void TSearch::DisplayPopulationStatistics(void)
{
	if (PopulationStatisticsDisplayFunction != NULL)
		(*PopulationStatisticsDisplayFunction)(Gen,BestPerf,AvgPerf,PerfVar);
	else {
		cout << "Generation " << Gen << ": Best = " << BestPerf;
		cout << ", Average = " << AvgPerf << ", Variance = " << PerfVar << endl;
	}
}


// Display the results of a search

void TSearch::DisplaySearchResults(void)
{
	if (SearchResultsDisplayFunction != NULL)
		(*SearchResultsDisplayFunction)(*this);
}


// Resume a search from a checkpoint file.  
// Note that this assumes that all function pointers in the current search object are the same
// as when the checkpoint file was saved, because function pointers are not saved in the 
// checkpoint file

void TSearch::ResumeSearch(void)
{
	// Restore the saved search object
    ReadCheckpointFile();
	// Restart the saved search
	DoSearch(1);
}


// Determine if the search is over

int TSearch::SearchTerminated(void)
{
	return (Gen >= MaxGens) || 
	       ((SearchTerminationFunction != NULL) && 
	        (*SearchTerminationFunction)(Gen,BestPerf,AvgPerf,PerfVar));
}


// **********
// Evaluation
// **********

// Evaluate the given vector
// Note that negative performances are treated as 0

double TSearch::EvaluateVector(TVector<double> &v, RandomState &rs)
{
	double perf = (*EvaluationFunction)(v, rs);
	
	return (perf<0)?0:perf;
}


// Evaluate a population range

void *EvaluatePopulationRange(void *arg)
{
  PopRangeSpec *prs = (PopRangeSpec *)arg;
  TSearch *s = prs->search;
  for (int i = prs->start; i <= prs->end; i++)
    s->Perf[i] = s->EvaluateVector(s->Population[i], s->RandomStates[i]);
  pthread_exit(NULL);
}


// Evaluate the current population, beginning with the STARTth individual

void TSearch::EvaluatePopulation(int start)
{
#ifdef THREADED_SEARCH  // Evaluate the population in parallel
  // Create threads
  if (THREAD_COUNT > 1) {
    int NumIndividuals = (PopulationSize() - start + 1)/THREAD_COUNT;
    pthread_t threads[THREAD_COUNT-1];
    PopRangeSpec psrs[THREAD_COUNT-1];
    int rc;
    for (int i = 1; i <= THREAD_COUNT - 1; i++) {
      psrs[i-1].search = this;
      psrs[i-1].start = (i-1)*NumIndividuals + start;
      psrs[i-1].end = i*NumIndividuals + start - 1;
      rc  = pthread_create(&threads[i-1], NULL, EvaluatePopulationRange, (void *)&psrs[i-1]);
      if (rc) {cerr << "Thread creation failed: " << rc << endl; exit(-1);}
    }
    // Evaluate the remaining individuals in the main thread
    for (int i = (THREAD_COUNT - 1)*NumIndividuals + start; i <= PopulationSize(); i++)
      Perf[i] = EvaluateVector(Population[i], RandomStates[i]);
    // Wait for all other threads to complete
    int status; 
    for (int i = 0; i <= THREAD_COUNT-2; i++)
      pthread_join(threads[i], (void **)&status);
  }
  else
    for (int i = start; i <= Population.Size(); i++) 
      Perf[i] = EvaluateVector(Population[i], RandomStates[i]);
  
#else // Evaluate the population serially
	for (int i = start; i <= Population.Size(); i++) 
		Perf[i] = EvaluateVector(Population[i], RandomStates[i]);
#endif
}


// *********
// Selection
// *********

// Compute the coefficients for linear fitness scaling.
// See Goldberg's book, pp. 76-79.

double LinearScaleFactor(double min, double max, double avg, double FMultiple)
{
	// Check that the scaled min will be greater than 0
	if (min > (FMultiple * avg - max)/(FMultiple - 1))
		// If so, do a full linear scaling
		{
			double delta = max - avg;
			if (delta > 0.0) return (FMultiple - 1) * avg/delta;
			else return 0.0;
		}
	else
		// Otherwise, scale as much as possible
		{
			double delta = avg - min;
			if (delta > 0.0) return avg/delta;
			else return 0.0;
		}
}


// Assign a normalized fitness to every individual in the population.  There are two methods:
// fitness proporationate and rank-based.  The rank-based method uses Baker's linear ranking 
// method (see Goldberg's book pp. 124-125 or Mitchell's book pp. 169-170).  The fitness 
// formula is derived as follows.  If the highest ranked individual (with rank 1) receives 
// MaxExpOffspring, then the fitness is given by y = m(x-1) + MaxExpOffspring.  Since the 
// sum of the fitness over all individuals must equal 1, we can apply this constraint to 
// solve for m in this linear equation.

void TSearch::UpdatePopulationFitness(void)
{
	int psize = PopulationSize();
	SortPopulation();
	switch (SelectMode) {
		// Calculate normalized fitness based on a fitness proportionate method with linear scaling
		case FITNESS_PROPORTIONATE:
			{
				double m = LinearScaleFactor(MinPerf,MaxPerf,AvgPerf,MaxExpOffspring);
				double total = 0;
				for (int i = 1; i <= psize; i++)
				{
					fitness[i] = m * (Perf[i] - AvgPerf) + AvgPerf;
					total = total + fitness[i];
				}
				for (int i = 1; i <= psize; i++)
					fitness[i] = fitness[i]/total;
				break;
			}
		// Calculate normalized fitness based on a rank-based method
		case RANK_BASED:
			for (int i = 1; i <= psize; i++)
				fitness[i] = (MaxExpOffspring + (2.0 - 2.0*MaxExpOffspring)*((i-1.0)/(psize-1)))/psize;
			break;
		default: cerr << "Invalid selection mode" << endl; exit(0);
	}
}


// *****************
// Genetic Operators
// *****************

// Gaussian mutation

void TSearch::MutateVector(TVector<double> &v)
{
	double magnitude;
	TVector<double> TempVector(1,vectorSize);

	// Generate a normally-distributed random magnitude
	magnitude = rs.GaussianRandom(0.0,MutationVar);
	// Generate a random unit vector
	rs.RandomUnitVector(TempVector);
	// Apply the mutation to V
	for (int i = 1; i <= vectorSize; i++)
		if (ConstraintVector[i])
			v[i] = clip(v[i] + magnitude * TempVector[i],MinSearchValue,MaxSearchValue);
		else
			v[i] = v[i] + magnitude * TempVector[i];
}


// Perform a modular uniform crossover between two individuals

void TSearch::UniformCrossover(TVector<double> &v1, TVector<double> &v2)
{
	if (crossPoints.Size() < 2) return;
	for (int i = 1; i <= crossPoints.Size() - 1; i++)
		if (ProbabilisticChoice(0.5)) 
			for (int j = crossPoints[i]; j < crossPoints[i+1]; j++) {
				double temp = v1[j];
				v1[j] = v2[j];
				v2[j] = temp;
			}
	if (ProbabilisticChoice(0.5)) 	
		for (int j = crossPoints[crossPoints.Size()]; j <= vectorSize; j++) {
			double temp = v1[j];
			v1[j] = v2[j];
			v2[j] = temp;
			}
}


// Perform a modular two-point crossover between two individuals

void TSearch::TwoPointCrossover(TVector<double> &v1, TVector<double> &v2)
{
	if (crossPoints.Size() < 2) return;
	int i1 = rs.UniformRandomInteger(1,crossPoints.Size());
	int i2 = i1;
	while (i2 == i1) 
		i2 = rs.UniformRandomInteger(1,crossPoints.Size());
	if (i1 > i2) {
		int t = i1;
		i1 = i2;
		i2 = t;
	}
	for (int i = crossPoints[i1]; i < crossPoints[i2]; i++) {
		double temp = v1[i];
		v1[i] = v2[i];
		v2[i] = temp;
	}
}


// ************
// Reproduction
// ************

// Create a new population

void TSearch::ReproducePopulation(void)
{
	switch (RepMode) {
		case HILL_CLIMBING: ReproducePopulationHillClimbing(); break;
		case GENETIC_ALGORITHM: ReproducePopulationGeneticAlgorithm(); break;
		default: cerr << "Invalid reproduction mode" << endl; exit(0);
	}
}


void TSearch::ReproducePopulationHillClimbing(void)
{
	int psize = PopulationSize();
	
	// Calculate population fitness
	UpdatePopulationFitness();	
	// Select the parents using Baker's stochastic universal sampling
	TVector<TVector<double> > ParentPopulation(1,psize);
	TVector<double> ParentPerf(1,psize);
	int j = 1;
	double sum = 0;
	double rand = rs.UniformRandom(0.0,1.0);
	for (int i = 1; (i <= psize) && (j <= psize); i++) {
		sum += psize * fitness[i];
		while (rand < sum) {
			ParentPopulation[j] = Population[i];
			ParentPerf[j] = Perf[i];
			j++;
			rand++;
		}
	}	
  // Replace the current population with the parent population
  Population = ParentPopulation;
	// If ReEvalFlag is set
	if (ReEvalFlag) {
    // reset BestPerf
    BestPerf = -1;
    // re-evaluate the parents
    EvaluatePopulation();
    // and update the performance values for the parents
    ParentPerf = Perf;
  }
  // Produce the new population by mutating each parent
  for (int i = 1; i <= psize; i++)
    MutateVector(Population[i]);
  // Evaluate the children
  EvaluatePopulation();
  // Restore each parent whose child's performance is worse
  for (int i = 1; i <= psize; i++)
    if (ParentPerf[i] > Perf[i]) {
      Population[i] = ParentPopulation[i];
      Perf[i] = ParentPerf[i];
    }
}


void TSearch::ReproducePopulationGeneticAlgorithm(void)
{
	int psize = PopulationSize();
	
	// Calculate population fitness
	UpdatePopulationFitness();
	// Determine the number of elite individuals in the new population
	int ElitePop = (int)floor(EFraction*psize + 0.5);
	// Select the rest of the population using Baker's stochastic universal sampling
	TVector<TVector<double> > TempPopulation = Population;
	int j = ElitePop+1;
	double sum = 0;
	double rand = rs.UniformRandom(0.0,1.0);
	for (int i = 1; (i <= psize) && (j <= psize); i++) {
		sum += (psize-ElitePop) * fitness[i];
		while (rand < sum) {
			Population[j++] = TempPopulation[i];
			rand++;
		}
	}	
	// Randomly shuffle the nonelite parents in preparation for crossover
  if (CrossProb > 0) {
    TVector<double> TempInd;
    for (int i = ElitePop+1; i <= psize; i++) {
      int k = rs.UniformRandomInteger(i,psize);
      TempInd = Population[k];
      Population[k] = Population[i];
      Population[i] = TempInd;
    }
  }
	// Apply mutation or crossover to each nonelite parent and compute the child's performance
	int i = ElitePop+1;
	TVector<double> Parent1, Parent2;
	while (i <= psize) {
		// Perform crossover with probability CrossProb
		if (ProbabilisticChoice(CrossProb) && (i < psize)) {
			Parent1 = Population[i];
			Parent2 = Population[i+1];
			switch (CrossMode) {
				case UNIFORM: UniformCrossover(Population[i],Parent2); break;
				case TWO_POINT: TwoPointCrossover(Population[i],Parent2); break;
				default: cerr << "Invalid crossover mode" << endl; exit(0);
			}
			// If the child is the same as the first parent after crossover, mutate it
			if (EqualVector(Population[i],Parent1)) MutateVector(Population[i]);
			i++;
		}
		// Otherwise, perform mutation
		else MutateVector(Population[i++]);
	}
  // Evaluate the new population
  if (ReEvalFlag) EvaluatePopulation();
  else EvaluatePopulation(ElitePop+1);
}


// Quicksort the population in descending order by performance

inline int partition(int first, int last, TVector<double> &perf, TVector<TVector<double> > &pop)
{
	int pivot = first;
	double pivot_value = perf[first];
	double temp1;
	TVector<double> temp2;
	
	for (int i = first; i <= last; i++) {
		if (perf[i] > pivot_value) {
			pivot++;
			if (i != pivot) {
				temp1 = perf[pivot]; perf[pivot] = perf[i]; perf[i] = temp1;
				temp2 = pop[pivot]; pop[pivot] = pop[i]; pop[i] = temp2;
			}
		}
	}
	temp1 = perf[pivot]; perf[pivot] = perf[first]; perf[first] = temp1;
	temp2 = pop[pivot]; pop[pivot] = pop[first]; pop[first] = temp2;
	
	return pivot;
}

inline void quicksort(int first, int last, TVector<double> &perf, TVector<TVector<double> > &pop)
{
	if (first < last) {
		int pivot = partition(first,last,perf,pop);
		quicksort(first,pivot-1,perf,pop);
		quicksort(pivot+1,last,perf,pop);
	}
}

void TSearch::SortPopulation(void)
{
	quicksort(1,Population.Size(),Perf,Population);
}


// ****************
// Input and Output
// ****************

// Read and write a search object
//
// Note that a complete representation of the state of a search cannot be stored
// in a file because of the function pointers. These i/o methods are primarily
// designed to support a simple checkpoint/restart facility.
//
// File format:
//  <Vector Size> <Population Size>
//  <Generation> <Max Generation>
//  <Random State>
//  <Selection Mode> <Reproduction Mode> <Crossover Mode>
//  <Search Initialized?> <Re-evaluation Flag> <Checkpoint Frequency>
//  <Search Constraint>
//  <Mutation Variance> <Crossover Prob>
//  <Crossover Template>
//  <Elitist Fraction> <Maximum Expected Offspring>
//  <Best Performance> <Best Vector>
//  <Performance 1> <Individual 1>
//  ...
//  <Performance N> <Individual N>
//  <RandomState 1>
//  ...
//  <RandomState N>

void TSearch::WriteCheckpointFile(void)
{
  ofstream bofs("search.cpt", ios::binary);
  int i;
  double d;
    
	// Write the vector size and population size
  bofs.write((const char*) &(vectorSize), sizeof(vectorSize));
  i = PopulationSize();
  bofs.write((const char*) &(i), sizeof(i));
	// Write the generation number and the maximum number of generations
  bofs.write((const char*) &(Gen), sizeof(Gen));
  bofs.write((const char*) &(MaxGens), sizeof(MaxGens));
	// Write the random state
	rs.BinaryWriteRandomState(bofs);                                           
	// Write the selection mode
	switch (SelectMode) {
		case FITNESS_PROPORTIONATE: i = 1; break;
		case RANK_BASED: i = 2; break;
		default: cerr << "Invalid selection mode" << endl; exit(0);
	}
  bofs.write((const char*) &(i), sizeof(i));      
	// Write the reproduction mode
	switch (RepMode) {
		case HILL_CLIMBING: i = 1; break;
		case GENETIC_ALGORITHM: i = 2; break;
		default: cerr << "Invalid reproduction mode" << endl; exit(0);
	}
  bofs.write((const char*) &(i), sizeof(i));
	// Write the crossover mode
	switch (CrossMode) {
		case UNIFORM: i = 1; break;
		case TWO_POINT: i = 2; break;
		default: cerr << "Invalid crossover mode" << endl; exit(0);
	}
  bofs.write((const char*) &(i), sizeof(i));  
	// Write the search initialized and re-evaluation flags, and the checkpoint frequency
  bofs.write((const char*) &(SearchInitialized), sizeof(SearchInitialized));
  bofs.write((const char*) &(ReEvalFlag), sizeof(ReEvalFlag));
  bofs.write((const char*) &(CheckpointInt), sizeof(CheckpointInt));
	// Write the search constraint vector
  ConstraintVector.BinaryWriteVector(bofs);
	// Write the mutation variance
  bofs.write((const char*) &(MutationVar), sizeof(MutationVar));
	// Write the crossover probability
  bofs.write((const char*) &(CrossProb), sizeof(CrossProb));
	// Write the crossover template
  crossTemplate.BinaryWriteVector(bofs);
	// Write the elitist fraction
	bofs.write((const char*) &(EFraction), sizeof(EFraction));
	// Write the max expected offspring
  bofs.write((const char*) &(MaxExpOffspring), sizeof(MaxExpOffspring));
	// Write out the peformance and parameter vector of the best individual
  bofs.write((const char*) &(BestPerf), sizeof(BestPerf));
  bestVector.BinaryWriteVector(bofs);
	// Write out the performance and parameter vector of each individual in the population
  for (int i = 1; i <= PopulationSize(); i++) {
    d = Performance(i);
    bofs.write((const char*) &(d), sizeof(d));
    Individual(i).BinaryWriteVector(bofs);
  }
  // Write out the random state for each individual in the population
  for (int i = 1; i <= PopulationSize(); i++)
    RandomStates[i].BinaryWriteRandomState(bofs);
}


void TSearch::ReadCheckpointFile(void)
{
  ifstream bifs("search.cpt", ios::binary);
  int i;
  double d;
  TVector<int> iv;
    
  // Read the vector size and population size
  bifs.read((char*) &(i), sizeof(i));
  SetVectorSize(i);
  bifs.read((char*) &(i), sizeof(i));
  SetPopulationSize(i);
	// Read the generation number and the maximum number of generations
  bifs.read((char*) &(Gen), sizeof(Gen));
  bifs.read((char*) &(i), sizeof(i));
  SetMaxGenerations(i);
	// Read the random state
	rs.BinaryReadRandomState(bifs);
	// Read the selection mode
	bifs.read((char*) &(i), sizeof(i));
	switch (i) {
    case 1: SetSelectionMode(FITNESS_PROPORTIONATE); break;
		case 2: SetSelectionMode(RANK_BASED); break;
		default: cerr << "Invalid selection mode" << endl; exit(0);
	}
	// Read the reproduction mode
	bifs.read((char*) &(i), sizeof(i));
	switch (i) {
		case 1: SetReproductionMode(HILL_CLIMBING);break;
		case 2: SetReproductionMode(GENETIC_ALGORITHM);break;
		default: cerr << "Invalid reproduction mode" << endl; exit(0);
	}
	// Read the crossover mode
  bifs.read((char*) &(i), sizeof(i));
	switch (i) {
		case 1: SetCrossoverMode(UNIFORM);break;
		case 2: SetCrossoverMode(TWO_POINT);break;
		default: cerr << "Invalid crossover mode" << endl; exit(0);
	}
	// Read the search initialized and re-evaluation flags, and the checkpoint frequency
  bifs.read((char*) &(SearchInitialized), sizeof(SearchInitialized));
  bifs.read((char*) &(ReEvalFlag), sizeof(ReEvalFlag));
  bifs.read((char*) &(CheckpointInt), sizeof(CheckpointInt));
  // Read the search constraint vector
  iv.BinaryReadVector(bifs);
  SetSearchConstraint(iv);
	// Read the mutation variance
  bifs.read((char*) &(d), sizeof(d));
  SetMutationVariance(d);
	// Read the crossover probability
  bifs.read((char*) &(d), sizeof(d));
	SetCrossoverProbability(d);
	// Read the crossover template
  iv.BinaryReadVector(bifs);
  SetCrossoverTemplate(iv);
	// Read the elitist fraciton
  bifs.read((char*) &(d), sizeof(d));
	SetElitistFraction(d);
	// Read the max expected offspring
  bifs.read((char*) &(d), sizeof(d));
	SetMaxExpectedOffspring(d);
	// Read the peformance and parameter vector of the best individual
  bifs.read((char*) &(BestPerf), sizeof(BestPerf));
  bestVector.BinaryReadVector(bifs);
	// Read the performance and parameter vector of each individual in the population
  for (int i = 1; i <= PopulationSize(); i++) {
    bifs.read((char*) &(d), sizeof(d));
    Perf[i] = d;
    Population[i].BinaryReadVector(bifs);
  }
  // Read in the random state for each individual in the populaton
  for (int i = 1; i <= PopulationSize(); i++)
    RandomStates[i].BinaryReadRandomState(bifs);
}
