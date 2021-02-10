//
//  StretchReceptor.cpp
//  one
//
//  Created by Eduardo on 9/26/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//

#include "StretchReceptor.h"

StretchReceptor::StretchReceptor(int nSegs, int nSR, double srvncgain, double srheadgain)
{
    SetStretchReceptorParams(nSegs, nSR, srvncgain, srheadgain);
}

void StretchReceptor::SetStretchReceptorParams(int nSegs, int nSR, double srvncgain, double srheadgain)
{
    NSEGS = nSegs;                  // Number of segments
    NSR = nSR;                      // Number of stretch receptors
    NSEGSSR = 6;                    // Number of segments that go into a stretch receptor
    SRvncgain = srvncgain;                // Stretch receptor gain
    SRheadgain = srheadgain;                // Stretch receptor gain

    NSEGSHEADSTART = 7;             // 7-12
    NSEGSHEAD = 14;                 // Number of segments for the sublateral head motorneurons
    NSEGSVNCSTART = 7;              // Segment where VNC starts

    normSegLenD.SetBounds(1, NSEGS);
    normSegLenV.SetBounds(1, NSEGS);
    D_sr.SetBounds(1, NSR);
    VA_sr.SetBounds(1, NSR);
    VP_sr.SetBounds(1, NSR);
}

void StretchReceptor::Update()
{
    double d, v;

    // Head
    d = 0.0;
    v = 0.0;
    for (int j = NSEGSHEADSTART; j < NSEGSHEADSTART + NSEGSHEAD; j++)
    {
        d += normSegLenD(j);
        v += normSegLenV(j);
    }
    HD_sr = SRheadgain*(d/NSEGSHEAD);
    HV_sr = SRheadgain*(v/NSEGSHEAD);

    // First four VC Neural Units (with three muscles each)
    for (int i = 1; i <= 6; i++){
        d = 0.0;
        v = 0.0;
        for (int j = 1; j <= NSEGSSR; j++)
        {
            d += normSegLenD(j+((i-1)*NSEGSSR)+NSEGSVNCSTART-1);
            v += normSegLenV(j+((i-1)*NSEGSSR)+NSEGSVNCSTART-1);
        }
        D_sr(i) = SRvncgain*(d/NSEGSSR);
        VA_sr(i) = SRvncgain*(v/NSEGSSR);

        d = 0.0;
        v = 0.0;
        for (int j = 1; j <= NSEGSSR; j++)
        {
            d += normSegLenD(j+((i-1)*NSEGSSR)+NSEGSVNCSTART-1+2);
            v += normSegLenV(j+((i-1)*NSEGSSR)+NSEGSVNCSTART-1+2);
        }
        VP_sr(i) = SRvncgain*(v/NSEGSSR);
    }
}
