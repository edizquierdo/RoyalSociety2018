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
    double chem_signal_scalar;
    int target_nrn_idx;
    bool enabled = false;

    void initialize(VecXY in_foodpos, double in_chem_signal_scalar, int in_target_nrn_idx)
    {
        foodpos = in_foodpos;
        chem_signal_scalar = in_chem_signal_scalar;
        target_nrn_idx = in_target_nrn_idx;
        enabled = true;
    }

    double get_CR_input(VecXY headpos)
    {
        // TODO: make this more accurate -- corner distance, or full diffusion sim
        return dist(headpos, foodpos) * chem_signal_scalar;
    }
};

