// ************************************************************
// A nervous system class (based on the CTRNN class)
//
// RDB
//  1/15 Created
// ************************************************************

#include "NervousSystem.h"
#include "random.h"
#include <stdlib.h>


// ****************************
// Constructors and Destructors
// ****************************

// The constructor

NervousSystem::NervousSystem(int newsize, int newmaxchemconns, int newmaxelecconns)
{
    SetCircuitSize(newsize, newmaxchemconns, newmaxelecconns);
}


// The destructor

NervousSystem::~NervousSystem()
{
    SetCircuitSize(0, 0, 0);
}


// *********
// Utilities
// *********

// Resize a circuit.

void NervousSystem::SetCircuitSize(int newsize, int newmaxchemconns, int newmaxelecconns)
{
    size = newsize;
    if (newmaxchemconns == -1) maxchemconns = size;
    else maxchemconns = min(newmaxchemconns, size);
    if (newmaxelecconns == -1) maxelecconns = maxchemconns;
    else maxelecconns = min(newmaxelecconns, maxchemconns);
    states.SetBounds(1,size);
    states.FillContents(0.0);
    paststates.SetBounds(1,size);  
    paststates.FillContents(0.0);
    outputs.SetBounds(1,size);
    outputs.FillContents(0.0);
    biases.SetBounds(1,size);
    biases.FillContents(0.0);
    gains.SetBounds(1,size);
    gains.FillContents(1.0);
    taus.SetBounds(1,size);
    taus.FillContents(1.0);
    Rtaus.SetBounds(1,size);
    Rtaus.FillContents(1.0);
    externalinputs.SetBounds(1,size);
    externalinputs.FillContents(0.0);
    NumChemicalConns.SetBounds(1,size);
    for (int i = 1; i <= size; i++)
        NumChemicalConns[i] = 0;
    chemicalweights.SetBounds(1,size,1,maxchemconns);
    NumElectricalConns.SetBounds(1,size);
    for (int i = 1; i <= size; i++)
        NumElectricalConns[i] = 0;
    electricalweights.SetBounds(1,size,1,maxelecconns);
    TempStates.SetBounds(1,size);
    TempOutputs.SetBounds(1,size);
    k1.SetBounds(1,size);
    k2.SetBounds(1,size);
    k3.SetBounds(1,size);
    k4.SetBounds(1,size);
}


// *********
// Accessors
// *********

double NervousSystem::ChemicalSynapseWeight(int from, int to)
{
    for (int i = 1; i <= NumChemicalConns(to); i++) {
        if (chemicalweights[to][i].from == from)
            return chemicalweights[to][i].weight;
    }
    return 0.0;
}


void NervousSystem::SetChemicalSynapseWeight(int from, int to, double value)
{
    // If the connection is already stored, just change its value
    for (int i = 1; i <= NumChemicalConns[to]; i++)
        if (chemicalweights[to][i].from == from) {
            chemicalweights[to][i].weight = value;
            return;
        };
    // Otherwise, make sure we have room for an additional connection ...
    if (NumChemicalConns[to] == maxchemconns) {
        cerr << "Maximum chemical synapses (" << maxchemconns << ") exceeded for neuron " << to << endl;
        exit(EXIT_FAILURE);
    }
    // ... and store it
    int i = ++NumChemicalConns[to];
    chemicalweights[to][i].from = from;
    chemicalweights[to][i].weight = value;
}


double NervousSystem::ElectricalSynapseWeight(int from, int to)
{
    for (int i = 1; i <= NumElectricalConns(to); i++) {
        if (electricalweights[to][i].from == from)
            return electricalweights[to][i].weight;
    }
    return 0.0;
}


void NervousSystem::InternalSetElectricalSynapseWeight(int from, int to, double value)
{
    // If the connection is already stored, just change its value
    for (int i = 1; i <= NumElectricalConns[to]; i++)
        if (electricalweights[to][i].from == from) {
            electricalweights[to][i].weight = value;
            return;
        };
    // Otherwise, make sure we have room for an additional connection ...
    if (NumElectricalConns[to] == maxelecconns) {
        cerr << "Maximum electrical synapses (" << maxelecconns << ") exceeded for neuron " << to << endl;
        exit(EXIT_FAILURE);
    }
    // ... and store it
    int i = ++NumElectricalConns[to];
    electricalweights[to][i].from = from;
    electricalweights[to][i].weight = value;
}

void NervousSystem::SetElectricalSynapseWeight(int n1, int n2, double value)
{
    if (value < 0) {
        cerr << "Electrical synapse weight between neurons " << n1 << " and " << n2 << " is negative: " << value << endl;
        exit(EXIT_FAILURE);
    }
    InternalSetElectricalSynapseWeight(n1, n2, value);
    InternalSetElectricalSynapseWeight(n2, n1, value);
}


// *******
// Control
// *******

// Randomize the states or outputs of a circuit.

void NervousSystem::RandomizeCircuitState(double lb, double ub)
{
    for (int i = 1; i <= size; i++)
        SetNeuronState(i, UniformRandom(lb, ub));
}

void NervousSystem::RandomizeCircuitState(double lb, double ub, RandomState &rs)
{
    for (int i = 1; i <= size; i++)
        SetNeuronState(i, rs.UniformRandom(lb, ub));
}

void NervousSystem::RandomizeCircuitOutput(double lb, double ub)
{
    for (int i = 1; i <= size; i++)
        SetNeuronOutput(i, UniformRandom(lb, ub));
}

void NervousSystem::RandomizeCircuitOutput(double lb, double ub, RandomState &rs)
{
    for (int i = 1; i <= size; i++)
        SetNeuronOutput(i, rs.UniformRandom(lb, ub));
}



// Integrate a circuit one step using Euler integration.
void NervousSystem::EulerStep(double stepsize)
{
    // Update past states (used for gap junctions)
    for (int i = 1; i <= size; i++){
        paststates[i] = states[i];
    }
    // Update the state of all neurons.
    for (int i = 1; i <= size; i++) {
        // External input
        double input = externalinputs[i];
        // Input from chemical synapses
        for (int j = 1; j <= NumChemicalConns[i]; j++)
            input += chemicalweights[i][j].weight * outputs[chemicalweights[i][j].from];
        // Input from electrical synapses
        for (int j = 1; j <= NumElectricalConns[i]; j++)
            input += electricalweights[i][j].weight * (paststates[electricalweights[i][j].from] - paststates[i]);
        // Take the step
        states[i] += stepsize * Rtaus[i] * (input - states[i]);
    }
    // Update the outputs of all neurons.
    for (int i = 1; i <= size; i++)
        outputs[i] = sigmoid(gains[i] * (states[i] + biases[i]));
}


// ****************
// Input and Output
// ****************

#include <iomanip>

ostream& operator<<(ostream& os, NervousSystem& c)
{
    // Set the precision
    os << setprecision(32);
    // Write the size, maxchemconns and maxelecconns
    os << c.size << " " << c.maxchemconns << " " << c.maxelecconns << endl << endl;
    // Write the time constants
    for (int i = 1; i <= c.size; i++)
        os << c.taus[i] << " ";
    os << endl << endl;
    // Write the biases
    for (int i = 1; i <= c.size; i++)
        os << c.biases[i] << " ";
    os << endl << endl;
    // Write the gains
    for (int i = 1; i <= c.size; i++)
        os << c.gains[i] << " ";
    os << endl << endl;
    // Write the chemical weights in sparse format (N from1 weight1 ... fromN weightN)
    for (int i = 1; i <= c.size; i++) {
        cout << c.NumChemicalConns[i] << "  ";
        for (int j = 1; j <= c.NumChemicalConns[i]; j++)
            os << c.chemicalweights[i][j].from << " " << c.chemicalweights[i][j].weight << "  ";
        os << endl;
    }
    os << endl;
    // Write the electrical weights in sparse format (N from1 weight1 ... fromN weightN)
    for (int i = 1; i <= c.size; i++) {
        cout << c.NumElectricalConns[i] << "  ";
        for (int j = 1; j <= c.NumElectricalConns[i]; j++)
            os << c.electricalweights[i][j].from << " " << c.electricalweights[i][j].weight << "  ";
        os << endl;
    }
    // Return the ostream
    return os;
}


istream& operator>>(istream& is, NervousSystem& c)
{
    // Read the sizes
    int size;
    is >> size;
    int maxchemconns;
    is >> maxchemconns;
    int maxelecconns;
    is >> maxelecconns;
    c.SetCircuitSize(size, maxchemconns, maxelecconns);
    // Read the time constants
    for (int i = 1; i <= size; i++) {
        is >> c.taus[i];
        c.Rtaus[i] = 1/c.taus[i];
    }
    // Read the biases
    for (int i = 1; i <= size; i++)
        is >> c.biases[i];
    // Read the gains
    for (int i = 1; i <= size; i++)
        is >> c.gains[i];
    // Read the chemical weights
    int n;
    for (int i = 1; i <= size; i++) {
        is >> n;
        for (int j = 1; j <= n; j++) {
            is >> c.chemicalweights[i][j].from;
            is >> c.chemicalweights[i][j].weight;
            c.NumChemicalConns[i]++;
        }
    }
    // Read the electrical weights
    for (int i = 1; i <= size; i++) {
        is >> n;
        for (int j = 1; j <= n; j++) {
            is >> c.electricalweights[i][j].from;
            is >> c.electricalweights[i][j].weight;
            c.NumElectricalConns[i]++;
        }
    }
    // Return the istream
    return is;
}
