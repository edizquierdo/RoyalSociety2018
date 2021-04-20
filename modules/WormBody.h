//
// WormBody.h
//
// An implmentation of the C. elegans body model described in the following paper (henceforth BBC):
//
// J.H. Boyle, S. Berri and N. Cohen (2012), Gait modulation in C. elegans:
// an integrated neuromechanical model, Front. Comput. Neurosci, 6:10
// doi: 10.3389/fncom.2012.00010
//
// Although this implementation is completely original, the body model itself should be identical to that of BBC
// with the exceptions commented as "NOTE n" within the code and explained at the top of the corresponding file.
//
//   NOTE 1: I fixed an obvious C macro bug in the original BBC code. In that code, NBAR is defined to be "NSEG+1",
//   but then used in, e.g., "2.0*NBAR", which expands into "2.0*NSEG+1" rather than "2.0*(NSEG+1)".
//
// Created by Randall Beer on 7/8/14.
// Copyright (c) 2014 Randall Beer. All rights reserved.
//
// Revision History

#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <cmath>

#include "Collide.h"

using namespace std;

// double M_PI = 3.1415926535897932384626433;

// If the symbol BBC_STRICT is defined, the deviations from the original BBC mentioned in NOTES 1 and 3
// are replaced with the orginals, so that the two models should be identical

//#define BBC_STRICT


// Settable constants

const double Medium             = 1.0;                           // Normalized medium drag coefficient (0 = water, 1 = agar)
const double L_worm             = 1.0e-3;                      // Length of worm in m
const int    N_segments         = 50;                    //YYY      // Number of segments
const double R_min              = 40.0e-6;                     // Minor radius of prolate ellipse body in m
const double C_agar_par_total   = 3.2e-3;                      // Total tangential drag coefficient for agar in kg/s
const double C_agar_perp_total  = 128e-3;                      // Total rod normal drag coefficient in agar in kg/s
const double C_water_par_total  = 3.3e-6;                      // Total rod tangential drag coefficient for water in kg/s
const double C_water_perp_total = 5.2e-6;                      // Total rod normal drag coefficient for water in kg/s
const double kappa_L            = (10.0e-3*N_segments)/24;     // Lateral spring constant in kg/s
const double kappa_D            = 350*kappa_L;                 // Diagonal spring constant in kg/s
const double kappa_M0           = 20*kappa_L;                  // Baseline active muscle spring constant in kg/s
const double beta_L             = 0.025*kappa_L;               // Lateral passive damping constant in s
const double beta_D             = 0.01*kappa_D;                // Diagonal passive damping constant in s
const double beta_M0            = 100*beta_L;                  // Baseline active damping constant in s
const double delta_M            = 0.65;                        // Rest muscle length scaling constant


// Derived constants

const int    N_rods       = N_segments+1;                            // Number of rods
const int    N_states     = 3*N_rods;                                // Total number of states in the body
const double L_seg        = L_worm/N_segments;                       // Length of an individual segment in m
const double D_min        = 2*R_min;                                 // Minor diameter of prolate ellipse body in m
#ifdef BBC_STRICT
const double C_agar_par   = C_agar_par_total/(2*N_segments + 1);     // Per rod tangential drag coefficient for agar in kg/s;  **** NOTE 1 ****
const double C_agar_perp  = C_agar_perp_total/(2*N_segments + 1);    // Per rod normal drag coefficient in agar in kg/s;       **** NOTE 1 ****
const double C_water_par  = C_water_par_total/(2*N_segments + 1);    // Per rod tangential drag coefficient for water in kg/s; **** NOTE 1 ****
const double C_water_perp = C_water_perp_total/(2*N_segments + 1);   // Per rod normal drag coefficient for water in kg/s;     **** NOTE 1 ****
#else
const double C_agar_par   = C_agar_par_total/(2*(N_segments + 1));   // Per rod tangential drag coefficient for agar in kg/s;  **** NOTE 1 ****
const double C_agar_perp  = C_agar_perp_total/(2*(N_segments + 1));  // Per rod normal drag coefficient in agar in kg/s;       **** NOTE 1 ****
const double C_water_par  = C_water_par_total/(2*(N_segments + 1));  // Per rod tangential drag coefficient for water in kg/s; **** NOTE 1 ****
const double C_water_perp = C_water_perp_total/(2*(N_segments + 1)); // Per rod normal drag coefficient for water in kg/s;     **** NOTE 1 ****
#endif
const double C_par        = (C_agar_par - C_water_par)*Medium + C_water_par;    // Per rod tangential drag coefficient in kg/s
const double C_perp       = (C_agar_perp - C_water_perp)*Medium + C_water_perp; // Per rod normal drag coefficient in kg/s




// Function prototypes

void InitializeBodyConstants(void);


// The WormBody class

class WormBody {
public:
    // Accessors
    inline double time() {return t;}
    inline double X(int i) {return Z[3*(i-1)];}     // YYY == SHOULD THIS BE i-1
    inline double Y(int i) {return Z[3*(i-1)+1];}   // YYY
    inline double Phi(int i) {return Z[3*(i-1)+2];} // YYY
    inline void SetDorsalSegmentActivation(int i, double a)
    {
        if (i > 0 && i <= N_segments)
            A_D_M[i-1] = fmax(0.0, fmin(a, 1.0));
        else {cerr << "SetDorsalSegmentActivation: " << i << " not in the range [1," << N_segments << "]" << endl; exit(EXIT_FAILURE);}
    }
    inline void SetVentralSegmentActivation(int i, double a)
    {
        if (i > 0 && i <= N_segments)
             A_V_M[i-1] = fmax(0.0, fmin(a, 1.0));
        else {cerr << "SetVentralSegmentActivation: " << i << " not in the range [1," << N_segments << "]" << endl; exit(EXIT_FAILURE);}
    }
    inline double DorsalSegmentLength(int i)
    {
        if (i > 0 && i <= N_segments) return L_D_L[i-1];
        else {cerr << "DorsalSegmentLength: " << i << " not in the range [1," << N_segments << "]" << endl; exit(EXIT_FAILURE);}
    }
    inline double VentralSegmentLength(int i)
    {
        if (i > 0 && i <= N_segments) return L_V_L[i-1];
        else {cerr << "VentralSegmentLength: " << i << " not in the range [1," << N_segments << "]" << endl; exit(EXIT_FAILURE);}
    }
    
    // YYY
    double RestingLength(int i);

    vector<CollisionObject> CollObjs;
    // void load_CollObjs(std::string collide_file) { CollObjs = load_objects(collide_file); }

    
    // Control
    void InitializeBodyState(double angle = 0.0, std::vector<CollisionObject> collObjs = std::vector<CollisionObject>());
    inline void StepBody(double h) {SemiImplicitBackwardEulerDAEStep(h);}

private:
    // More control
    inline void F(int start = 0, int end = N_rods)
    {
        UpdateKinematics(start,end);
        UpdateForces(start, end);
        UpdateResiduals(start,end);
    }
    void UpdateKinematics(int start = 0, int end = N_rods);
    void UpdateForces(int start = 0, int end = N_rods);
    void UpdateResiduals(int start = 0, int end = N_rods);
    void SemiImplicitBackwardEulerDAEStep(double h);
    void NumericaldFdZ(double J[N_states][N_states]);
    void NumericaldFdZp(double J[N_states][N_states]);
    void LinearSolve(double M[N_states][N_states], double dZ[], double deltaZ[]);
    
    // Debugging
    void DebugPrintVector(const char s[], double *v, int n = N_states)
    {
        cout << s;
        for (int i = 0; i < n; i++) cout << v[i] << " ";
        cout << endl;
    };
    void DebugPrintMatrix(const char s[], double M[N_states][N_states])
    {
        cout << s;
        for (int i = 0; i < N_states; i++) {
            for (int j = 0; j < N_states; j++) {
                cout << M[i][j] << " ";
            }
            cout << endl;
        }
    };

    
    
    // Instance variables
    double t;
    double Z[N_states], dZ[N_states], Residuals[N_states];
    double A_D_M[N_segments],A_V_M[N_segments];
    double sinPhi[N_rods], cosPhi[N_rods];
    double L_D_D[N_segments],dL_D_D[N_segments],L_V_D[N_segments],dL_V_D[N_segments];
    double uD_D_x[N_segments],uD_D_y[N_segments],uD_V_x[N_segments],uD_V_y[N_segments];
    double L_D_L[N_segments],dL_D_L[N_segments], L_V_L[N_segments], dL_V_L[N_segments];
    double uL_D_x[N_segments],uL_D_y[N_segments],uL_V_x[N_segments],uL_V_y[N_segments];
    double f_D_x[N_rods],f_D_y[N_rods],f_V_x[N_rods],f_V_y[N_rods];    
};

