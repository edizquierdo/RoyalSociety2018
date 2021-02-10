//
// LocomotionCircuit.h
//
//  Copyright (c) 2015 Eduardo Izquierdo. All rights reserved.
//

#include "VectorMatrix.h"
#include "random.h"
#include <cmath>
using namespace std;

class Muscles {
public:

    Muscles(int nmuscles = 24, double t_muscle = 0.1);

    void SetMuscleParams(int nmuscles, double t_muscle);

    void InitializeMuscleState();

    void EulerStep(double StepSize);

    void SetDorsalMuscleInput(int muscle, double input){V_input[muscle][1] = input;};
    void SetVentralMuscleInput(int muscle, double input){V_input[muscle][2] = input;};

    void SetDorsalMuscleActivation(int muscle, double activation){V_muscle[muscle][1] = activation;};
    void SetVentralMuscleActivation(int muscle, double activation){V_muscle[muscle][2] = activation;};

    double DorsalMuscleOutput(int muscle){return V_muscle[muscle][1];};
    double VentralMuscleOutput(int muscle){return V_muscle[muscle][2];};

    TMatrix<double> V_muscle;
    TMatrix<double> V_input;

    double T_muscle;
    int Nmuscles;
};
