//
//  Muscles.cpp
//  one
//
//  Created by Eduardo Izquierdo on 3/26/15.
//  Copyright (c) 2015 Eduardo Izquierdo. All rights reserved.
//

#include "Muscles.h"

Muscles::Muscles(int nmuscles, double t_muscle)
{
    SetMuscleParams(nmuscles, t_muscle);
}

void Muscles::SetMuscleParams(int nmuscles, double t_muscle)
{
    Nmuscles = nmuscles;
    T_muscle = t_muscle;
    V_muscle.SetBounds(1,Nmuscles,1,2);
    V_input.SetBounds(1,Nmuscles,1,2);
}

void Muscles::InitializeMuscleState()
{
    V_muscle.FillContents(0.0);
    V_input.FillContents(0.0);
}

void Muscles::EulerStep(double StepSize)
{
    for (int i = 1; i <= Nmuscles; i++)
    {
        V_muscle[i][1] += StepSize*((V_input[i][1] - V_muscle[i][1])/T_muscle);
        V_muscle[i][2] += StepSize*((V_input[i][2] - V_muscle[i][2])/T_muscle);
    }
}