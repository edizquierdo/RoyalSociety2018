// *******************************
// Various Random Number Utilities
//
// RDB 2/95
// *******************************

#pragma once

#include "VectorMatrix.h"
#include <fstream>

using namespace std;

// Global defines for ran1

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


// Functions to manipulate the global random state for backward compatibility

void SetRandomSeed(long seed);
long GetRandomSeed(void);
void WriteRandomState(ostream& os);
void BinaryWriteRandomState(ofstream& bofs);
void ReadRandomState(istream& is);
void BinaryReadRandomState(ifstream& bifs);
double UniformRandom(double min,double max);
int UniformRandomInteger(int min,int max);
double GaussianRandom(double mean, double variance);
void RandomUnitVector(TVector<double> &v);
int ProbabilisticChoice(double prob);


// The RandomState class declaration

class RandomState {
public:
  // The constructor
  RandomState(long seed = 0) {SetRandomSeed(seed); gaussian_flag = 0;};
  // The destructor
  ~RandomState() {};
  
  // Accessors
  void SetRandomSeed(long seed);
  long GetRandomSeed(void);
  
  // Helper functions
  double ran1(void);
  void GenerateNormals(void);
  
  // Return random deviates
  double UniformRandom(double min,double max);
  int UniformRandomInteger(int min,int max);
  double GaussianRandom(double mean, double variance);
  void RandomUnitVector(TVector<double> &v);
  int ProbabilisticChoice(double prob);
  
  // Input/Output 
  void WriteRandomState(ostream& os);
  void BinaryWriteRandomState(ofstream& bofs);
  void ReadRandomState(istream& is);
  void BinaryReadRandomState(ifstream& bifs);
  

  long seed, idum, iy, iv[NTAB];
  int gaussian_flag;
  double gX1, gX2;
};
