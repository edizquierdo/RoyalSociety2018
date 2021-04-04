// ************************************************************
// A nervous system class (based on the CTRNN class)
//
// RDB 
//  1/15 Created
// ************************************************************

#include "VectorMatrix.h"
#include "random.h"

#include "packages/json.hpp"

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <algorithm>

#pragma once

// #define CONNFILE "data/connectome.csv"
// #define NAMESMAP "data/names_map.csv"

using json = nlohmann::json;

// An entry in a sparse weight matrix

struct weightentry {int from; double weight;};


// The sigmoid function

inline double sigmoid(double x)
{
  return 1/(1 + exp(-x));
}


// The inverse sigmoid function

inline double InverseSigmoid(double y)
{
  return log(y/(1-y));
}

// utility functions for computing size and maximum connections
int compute_size(json & neurons)
{
    return std::distance(neurons.begin(), neurons.end());
}

int compute_maxconn(json & connections, string conn_type)
{
    std::unordered_map<string,int> counts = std::unordered_map<string,int>();

    // count up for for each connection
    for (auto & conn : connections)
    {
        if (conn["type"] == conn_type)
        {
            // if it is of the proper type, iterate the counter
            // REVIEW: target or presynaptic neuron?
            auto it = counts.find(conn["to"]);
            if (it != counts.end())
            {
                it->second += 1;
            }
            else 
            {
                // create the counter if it does not exist
                counts[conn["to"]] = 1;
            }
        }
    }

    // return max element of counts array
    return std::max_element(counts.begin(), counts.end())->second;
}



// The NervousSystem class declaration

class NervousSystem {
    public:
        // The constructor
        NervousSystem(int size = 0, int maxchemconns = -1, int maxelecconns = -1);
        // json based ctor
        NervousSystem(json & ns_data);
        // The destructor
        ~NervousSystem();
        
        // Accessors
        int CircuitSize(void) {return size;};
        void SetCircuitSize(int newsize, int maxchemconns, int maxelecconns);
        double NeuronState(int i) {return states[i];};
        void SetNeuronState(int i, double value) {states[i] = value;outputs[i] = sigmoid(gains[i]*(states[i] + biases[i]));};
        double NeuronOutput(int i) {return outputs[i];};
        void SetNeuronOutput(int i, double value) {outputs[i] = value; states[i] = InverseSigmoid(value)/gains[i] - biases[i];};
        double NeuronBias(int i) {return biases[i];};
        void SetNeuronBias(int i, double value) {biases[i] = value;};
        double NeuronGain(int i) {return gains[i];};
        void SetNeuronGain(int i, double value) {gains[i] = value;};
        double NeuronTimeConstant(int i) {return taus[i];};
        void SetNeuronTimeConstant(int i, double value) {taus[i] = value; Rtaus[i] = 1/value;};
        double NeuronExternalInput(int i) {return externalinputs[i];};
        void SetNeuronExternalInput(int i, double value) {externalinputs[i] = value;};
        double ChemicalSynapseWeight(int from, int to);
        void SetChemicalSynapseWeight(int from, int to, double value);
        double ElectricalSynapseWeight(int from, int to);
        void InternalSetElectricalSynapseWeight(int from, int to, double value);
        void SetElectricalSynapseWeight(int n1, int n2, double value);

        // Input and output
        friend ostream& operator<<(ostream& os, NervousSystem& c);
        friend istream& operator>>(istream& is, NervousSystem& c);
                            
        // Control
        void RandomizeCircuitState(double lb, double ub);
        void RandomizeCircuitState(double lb, double ub, RandomState &rs);
        void RandomizeCircuitOutput(double lb, double ub);
        void RandomizeCircuitOutput(double lb, double ub, RandomState &rs);
        void EulerStep(double stepsize);
        //void RK4Step(double stepsize);

        // loading
        // TODO: deprecate
        std::unordered_map<string,int> namesMap;
        void load_connectome(string connfile);
        void load_namesMap(string namesMap_file);
        
        void loadJSON_neurons(json & neurons);
        void loadJSON_connections(json & conn);
		
        int size, maxchemconns, maxelecconns;
        TVector<double> states, outputs, biases, gains, taus, Rtaus, externalinputs;
        TVector<double> paststates;  
        TVector<int> NumChemicalConns, NumElectricalConns;
        TMatrix<weightentry> chemicalweights, electricalweights;
        TVector<double> TempStates,TempOutputs,k1,k2,k3,k4;
        std::unordered_map<string,double> wgts_NMJ;
};

