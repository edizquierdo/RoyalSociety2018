// *******************************
// Various Random Number Utilities
//
// RDB 2/95
// *******************************

#include "random.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;


// A global random state

RandomState GRS;

// Functions to manipulate the global random state for backward compatibility

void SetRandomSeed(long seed) {GRS.SetRandomSeed(seed);};
long GetRandomSeed(void) {return GRS.GetRandomSeed();};
void WriteRandomState(ostream& os) {GRS.WriteRandomState(os);};
void BinaryWriteRandomState(ofstream& bofs) {GRS.BinaryWriteRandomState(bofs);};
void ReadRandomState(istream& is) {GRS.ReadRandomState(is);};
void BinaryReadRandomState(ifstream& bifs) {GRS.BinaryReadRandomState(bifs);};
double UniformRandom(double min,double max) { return GRS.UniformRandom(min, max);};
int UniformRandomInteger(int min,int max) {return GRS.UniformRandomInteger(min, max);};
double GaussianRandom(double mean, double variance) {return GRS.GaussianRandom(mean, variance);};
void RandomUnitVector(TVector<double> &v) {GRS.RandomUnitVector(v);};
int ProbabilisticChoice(double prob) {return GRS.ProbabilisticChoice(prob);};


// Return a uniform deviate between 0 and 1 exclusive with a period >10^8.
// Adapted from Ran1 in Numerical Recipes, p. 280.

double RandomState::ran1(void)
{
	if (!iy) SetRandomSeed(1);
	long k = idum/IQ;
	idum = IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	int j = iy/NDIV;   
	iy = iv[j];
	iv[j] = idum;
	double temp = AM * iy;
	if (temp > RNMX) return RNMX;
	else return temp;
}


// Seed the random number generator.
// Adapted from Ran1 from Numerical Recipes, p. 280.

void RandomState::SetRandomSeed(long s)
{
	idum = seed = s;
	if (idum == 0) idum = 1;
	for (int j = NTAB+7; j >= 0; j--) {
		long k = (idum)/IQ;
		idum = IA*(idum-k*IQ)-IR*k;
		if (idum < 0) idum += IM;
		if (j < NTAB) iv[j] = idum;
	}
	iy = iv[0];
}


// Return the current random seed

long RandomState::GetRandomSeed(void)
{
	return seed;
}


// Write a random state to a stream

void RandomState::BinaryWriteRandomState (ofstream& bosf)
{
  bosf.write((const char*) &(seed), sizeof(seed));
  bosf.write((const char*) &(idum), sizeof(idum));
  bosf.write((const char*) &(iy), sizeof(iy));
  bosf.write((const char*) &(gaussian_flag), sizeof(gaussian_flag));
  bosf.write((const char*) &(gX1), sizeof(gX1));
  bosf.write((const char*) &(gX2), sizeof(gX2));
  for (int i = 0; i < NTAB; i++)
    bosf.write((const char*) &(iv[i]), sizeof(iv[i]));
}

void RandomState::WriteRandomState(ostream& os)
{
	os << seed << " " << idum << " " << iy << " ";
  os << gaussian_flag << " " << gX1 << " " << gX2 << " ";
	for (int i = 0; i < NTAB; i++)
		os << iv[i] << " ";
	os << endl;
}


// Read a random state from a stream

void RandomState::ReadRandomState(istream& is)
{
	is >> seed;
	is >> idum;
	is >> iy;
  is >> gaussian_flag;
  is >> gX1;
  is >> gX2;
	for (int i = 0; i < NTAB; i++)
		is >> iv[i];
}

void RandomState::BinaryReadRandomState (ifstream& bisf)
{
  bisf.read((char*) &(seed), sizeof(seed));
  bisf.read((char*) &(idum), sizeof(idum));
  bisf.read((char*) &(iy), sizeof(iy));
  bisf.read((char*) &(gaussian_flag), sizeof(gaussian_flag));
  bisf.read((char*) &(gX1), sizeof(gX1));
  bisf.read((char*) &(gX2), sizeof(gX2));
  for (int i = 0; i < NTAB; i++)
    bisf.read((char*) &(iv[i]), sizeof(iv[i]));
}


// Return a uniformly-distributed random double between MIN and MAX exclusive.

double RandomState::UniformRandom(double min, double max)
{
	return (max - min) * ran1() + min;
}


// Return a uniformly-distributed random integer between MIN and MAX inclusive.

int RandomState::UniformRandomInteger(int min, int max)
{
	return (int)floor(0.5+UniformRandom(min-0.5,max+0.5));
}


// Generate two normally-distributed random variables gx1 and gX2 for use by
// GaussianRandom.  Based on the algorithm described in Volume 2 of "The Art
// of Computer Progamming" by Donald Knuth (p. 117).

void RandomState::GenerateNormals(void)
{
	double v1,v2,s,d;

	do
		{
			v1 = UniformRandom(-1.0,1.0);
			v2 = UniformRandom(-1.0,1.0);
			s = v1 * v1 + v2 * v2;
		}
	while (s >= 1.0 || s == 0.0);

	d = sqrt((-2 * log(s))/s);

	gX1 = v1 * d;
	gX2 = v2 * d;
}


// Generate a Gaussian random variable.

double RandomState::GaussianRandom(double mean, double variance)
{
	if (!gaussian_flag) {
		GenerateNormals();
		gaussian_flag = 1;
		return(sqrt(variance) * gX1 + mean);
	}
	else {
		gaussian_flag = 0;
		return(sqrt(variance) * gX2 + mean);
	}
}


// Generate a random unit vector.  This works by first generating a vector
// each of whose elements is a random Gaussian and then normalizing the
// resulting vector.  See Volume 2 of "The Art of Computer Programming"
// by Donald Knuth (pp. 130-131).

void RandomState::RandomUnitVector(TVector<double> &v)
{
	double r = 0.0;

	for (int i = v.LowerBound(); i <= v.UpperBound(); i++)
	{
		v[i] = GaussianRandom(0,1);
		r += v[i] * v[i];
	}
	r = sqrt(r);
	for (int i = v.LowerBound(); i <= v.UpperBound(); i++)
		v[i] = v[i] / r;
}


// Return 1 with a probability PROB, else return 0

int RandomState::ProbabilisticChoice(double prob)
{
	return (UniformRandom(0.0,1.0) <= prob)?1:0;
}
