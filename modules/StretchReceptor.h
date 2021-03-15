//
//  StretchReceptor.hpp
//  one
//
//  Created by Eduardo on 9/26/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//

#include "VectorMatrix.h"
#include "random.h"
#include <cmath>
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

