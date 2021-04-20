//
//  StretchReceptor.hpp
//  one
//
//  Created by Eduardo on 9/26/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//

#include <cmath>

#include "VectorMatrix.h"
#include "random.h"
#include "Collide.h"

using namespace std;

class StretchReceptor {
public:
    
    StretchReceptor(int nSegs = 50, int nSR = 7, double SRgain = 0.0, double SRHeadgain = 0.0);

    void SetStretchReceptorParams(int nSegs, int nSR, double SRVNCgain, double SRHeadgain);
    
    void Update();
    
    void SetDorsalInput(int seg, double normlen){normSegLenD(seg) = normlen;};
    void SetVentralInput(int seg, double normlen){normSegLenV(seg) = normlen;};
    

    double HeadDorsalOutput(){return HD_sr;};
    double HeadVentralOutput(){return HV_sr;};
    double VCDorsalOutput(int i){return D_sr(i);};
    double VCVentralAOutput(int i){return VA_sr(i);};
    double VCVentralPOutput(int i){return VP_sr(i);};
    
    double NSR;
    double NSEGS;
    double NSEGSSR;
    double SRvncgain;
    double SRheadgain;
    double NSEGSHEADSTART,NSEGSHEAD, NSEGSVNCSTART;
    TVector<double> normSegLenD;
    TVector<double> normSegLenV;
    double HD_sr;
    double HV_sr;
    TVector<double> D_sr;
    TVector<double> VA_sr;
    TVector<double> VP_sr;
    
};


class ChemoReceptor {
public:
    VecXY foodpos;
    int target_nrn_idx;
    bool enabled = false;
    // AWA consts
    double alpha; double beta; double gamma;

    // storing value of fast and slow sense
    double F_i; double S_i;
    double F_im1; double S_im1;

    double C_ixy; double out_AWA_stim;

    void initialize(VecXY in_foodpos, int in_target_nrn_idx, double in_alpha, double in_beta, double in_gamma)
    {
        enabled = true;
        foodpos = in_foodpos;
        target_nrn_idx = in_target_nrn_idx;

        alpha = in_alpha; 
        beta = in_beta; 
        gamma = in_gamma;

        F_im1 = 0.0;
        S_im1 = 0.0;
    }

    double get_concentration(VecXY headpos)
    {
        // TODO: make this more accurate -- corner distance, or full diffusion sim
        return dist(headpos, foodpos);
    }


    /* 
    same function as `comp_sensory()` from 
    https://github.com/mivanit/CE_learning/blob/main/CE_learn/wormSim.py
    */
    double comp_sensory(VecXY headpos, double StepSize)
    {
        C_ixy = get_concentration(headpos);

        // iterate fast and slow sense
		F_i = F_im1 + StepSize * ( (alpha * C_ixy) - (beta * F_im1) );
		S_i = S_im1 + StepSize * ( gamma * (F_im1 - S_im1) );

        // update stim for AWA neuron
        out_AWA_stim = F_i - S_i;
        
        // update prev-timestep values
        F_im1 = F_i; S_im1 = S_i;

        return out_AWA_stim;
    }
};

