//
//  Worm.hpp
//  one
//
//  Created by Eduardo Izquierdo on 9/25/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//

#include "VectorMatrix.h"
#include "random.h"
#include "WormBody.h"
#include "NervousSystem.h"
#include "Muscles.h"
#include "StretchReceptor.h"
#include "util.h"

#include <cmath>
#include <string>

#define PI 3.14159265

using namespace std;

// Parameters
const int N_muscles = 24;           // Number of muscles alongside the body
const int N_units = 6;              // Number of neural units
const int N_neuronsperunit = 6;     // Number of neurons in a neural unit

const int N_stretchrec = 6;         // N_units + 1 // Number of stretch receptors
//
const double T_muscle = 0.1;        // Muscle time constant

const int HeadMotorNeuronMuscles = 6;  // Head motorneurons innervate first 8 muscles (temporarily first 6)
const int VNCMuscleStart = 7;           // VNC motorneurons innervate starting from 7th muscle
const int NmusclePerNU = 3;             // All the way down to 24, in groups of 3 per unit


// Body segment name conventions
const int Head = 1;
const int Tail = N_segments;

class Worm {
public:
    
    Worm(TVector<double> &v, double output);
    Worm(json & params);

    void InitializeState(RandomState &rs, double angle, std::string collide_file);
    void HeadStep(double StepSize, double output);
    void Step(double StepSize, double output); 
    
    void DumpBodyState(ofstream &ofs, int skips);
    void DumpActState(ofstream &ofs, int skips);
    void DumpActState_header(ofstream &ofs);
    void DumpVoltage(ofstream &ofs, int skips);
    void DumpParams(ofstream &ofs);
    
    double CoMx();
    double CoMy();
    void Curvature(TVector<double> &c);
    double Orientation();
    
    WormBody b;
    Muscles m;
    NervousSystem n;
    StretchReceptor sr;
    ChemoReceptor chemo_re;
    NervousSystem h;
    
    double t; // Time
    
    // Neuromuscular junctions
    double NMJ_DB, NMJ_VBa, NMJ_VBp, NMJ_DD, NMJ_VDa, NMJ_VDp;
    double NMJ_SMDD, NMJ_RMDD, NMJ_SMDV, NMJ_RMDV;
    double NMJ_Gain_Map;

    // Neuron name conventions
    int DB = -1;
    int DD = -1;
    int VBA = -1;
    int VDA = -1;
    int VBP = -1;
    int VDP = -1;

    int SMDD = -1;
    int RMDD = -1;
    int SMDV = -1;
    int RMDV = -1;
    
    TVector<double> NMJ_Gain;
    
    // Head oscillator
    double dorsalinput1, ventralinput1, dorsalinput2, ventralinput2;
    double headFreq, headDelay, headGain, headBias;
};
