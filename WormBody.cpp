//
// WormBody.cpp
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
//   NOTE 2: There is a discrepency between the BBC code and paper as to whether or not the factor of 2 here
//   should be raised to the fourth power (16). I have followed the BBC code in raising it to the fourth power
//   (Boyle, personal communication, June 30, 2014).
//
//   NOTE 3: There is a major discrepency between the BBC code on the one hand and the BBC paper and Boyle thesis
//   on the other as to the formula for dphi/dt. The BBC code uses (f_odd_par/C_par)/(M_PI*2.0*R[i]), but both documents
//   use 2*f_odd_par/(R[i]*C_par). The argument for the latter formula makes perfect sense to me (see Boyle Thesis,
//   p. 76, 128). However, I cannot understand where the first forumula comes from and Boyle was unable to shed any
//   light on the matter (Boyle, personal communication, June 30, 2014). Which formula used does make a small but nontrivial
//   difference to the results. I have chosen to use the formula from the publications.
//
// Created by Randall Beer on 7/8/14.
// Copyright (c) 2014 Randall Beer. All rights reserved.
//
// Revision History
//  6/30/15 - Fixed bounds errors in UpdateForces (thanks to Eduardo Izquierdo)

#include "WormBody.h"
#include <cfloat>

using namespace std;


// Global constants

double R[N_rods];                     // Rod radii in m
double L_D0[N_rods];                  // Rest length of each diagonal element in m
double L_L0[N_rods];                  // Rest length of each lateral element in m
double L_min[N_rods];                 // Minimal length of each lateral muscle in m
double L_L0_minus_L_min[N_rods];      // Precomputed difference between the above two vectors


// Initialize various vectors of per rod body constants
// This only needs to be done ONCE, not once per instance

void InitializeBodyConstants(void)
{
    double NS2 = N_segments/2.0, LS2 = L_seg*L_seg, r;
    for (int i = 0; i < N_rods; i++)
        R[i] = R_min * fabs(sin(acos((i - NS2)/(NS2 + 0.2))));
    for (int i = 0; i < N_segments; i++) {
        r = R[i] - R[i+1]; L_L0[i] = sqrt(LS2+r*r);
        r = R[i] + R[i+1]; L_D0[i] = sqrt(LS2+r*r);
        L_min[i] = L_L0[i] * (1 - delta_M*r/D_min);
        L_L0_minus_L_min[i] = L_L0[i] - L_min[i];
    }
}


// Initialize the state of the body
// Note that, for a DAE, initial states cannot be chosen arbitrarily. They must be consistent with the DAE.
// Note also that the integration method SemiImplicitBackwardEulerDAEStep assumes that the initial dZ is 0

void WormBody::InitializeBodyState(void)
{
    t = 0.0;
    for (int i = 0; i < N_rods; i++) {
        int i3 = 3*i;
        Z[i3] = i*L_seg; Z[i3+1] = 0.0; Z[i3+2] = M_PI/2;
        dZ[i3] = dZ[i3+1] = dZ[i3+2] = 0.0;
    }
    UpdateKinematics();
    for (int i = 0; i < N_segments; i++)
        A_D_M[i] = A_V_M[i] = 0.0;
}

// YYY
double WormBody::RestingLength(int i)
{
    return L_L0[i-1];
}

// Update the body states one step using the semi-implicit backward Euler method -[(1/h)*dF/dZ' + dF/dZ].(Z_n+1 - Z_n) = F(Z_n,0)
// Essentially, this is backward Euler or 1st-order BDF, but with only one Newton iteration per step
// When F(Z,dZ) = 0 can be written as dZ = f(Z), this reduces to [I/h - df/dZ].(Z_n+1 - Z_n) = f(Z_n)

void WormBody::SemiImplicitBackwardEulerDAEStep(double h)
{
    double deltaZ[N_states], dFdZ[N_states][N_states] = {{0}}, M[N_states][N_states] = {{0}};
    
    // Update the residuals (assumes kinematics is up to date)
    UpdateForces();
    UpdateResiduals();
    // Compute and save dF/dZ
    NumericaldFdZ(dFdZ);
    // Compute and save dF/dZ'
    NumericaldFdZp(M);
    // Form the matrix M = -[(1/h)*dF/dZ' + dF/dZ]
    double h1 = 1.0/h;
    for (int i = 0; i < N_states; i++)
        for (int j = 0; j < N_states; j++)
            M[i][j] = -(h1*M[i][j] + dFdZ[i][j]);
    // Solve M . deltaZ = F for deltaZ
    LinearSolve(M, Residuals, deltaZ);
    // Compute Z_n+1 = Z_n + deltaZ
    for (int i = 0; i < N_states; i++) Z[i] += deltaZ[i];
    // Bring the kinematics up to date
    UpdateKinematics();
    // Increment time
    t += h;
}


// Update the kinematics of the body given the current state

void WormBody::UpdateKinematics(int start, int end)
{
    int start2 = max(0,start-1), end2 = min(end+1,N_rods);
    double p_D_x[N_rods],p_D_y[N_rods],p_V_x[N_rods],p_V_y[N_rods];
    double dp_D_x[N_rods],dp_D_y[N_rods],dp_V_x[N_rods],dp_V_y[N_rods];

    // Update rod endpoint positions and velocities
    for (int i = start2; i < end2; i++) {
        int i3 = 3*i;
        
        double x = Z[i3], y = Z[i3+1], phi = Z[i3+2];
        sinPhi[i] = sin(phi);
        cosPhi[i] = cos(phi);
        p_D_x[i] = x + R[i]*cosPhi[i];
        p_D_y[i] = y + R[i]*sinPhi[i];
        p_V_x[i] = x - R[i]*cosPhi[i];
        p_V_y[i] = y - R[i]*sinPhi[i];
        
        double dx = dZ[i3], dy = dZ[i3+1], dphi = dZ[i3+2];
        double
        arm = R[i]*dphi,
        dxtemp = -arm*sinPhi[i],  // cos(Phi + PI/2) <=> -sin(Phi)
        dytemp = arm*cosPhi[i];   // sin(Phi + PI/2) <=>  cos(Phi)
        dp_D_x[i] = dx + dxtemp;
        dp_D_y[i] = dy + dytemp;
        dp_V_x[i] = dx - dxtemp;
        dp_V_y[i] = dy - dytemp;
    }
    
    // Update element lengths, directions and velocities
    for (int i = start2; i < end2 - 1; i++) {
        int i1 = i+1;
        double D_x, D_y, V_x, V_y;
        // Diagonal elements
        D_x = p_V_x[i1] - p_D_x[i];
        D_y = p_V_y[i1] - p_D_y[i];
        V_x = p_D_x[i1] - p_V_x[i];
        V_y = p_D_y[i1] - p_V_y[i];
        L_D_D[i] = sqrt(D_x*D_x + D_y*D_y);
        uD_D_x[i] = D_x/L_D_D[i];
        uD_D_y[i] = D_y/L_D_D[i];
        dL_D_D[i] = (dp_V_x[i1] - dp_D_x[i])*uD_D_x[i] + (dp_V_y[i1] - dp_D_y[i])*uD_D_y[i];
        L_V_D[i] = sqrt(V_x*V_x + V_y*V_y);
        uD_V_x[i] = V_x/L_V_D[i];
        uD_V_y[i] = V_y/L_V_D[i];
        dL_V_D[i] = (dp_D_x[i1] - dp_V_x[i])*uD_V_x[i] + (dp_D_y[i1] - dp_V_y[i])*uD_V_y[i];
        // Lateral elements
        D_x = p_D_x[i1] - p_D_x[i];
        D_y = p_D_y[i1] - p_D_y[i];
        V_x = p_V_x[i1] - p_V_x[i];
        V_y = p_V_y[i1] - p_V_y[i];
        L_D_L[i] = sqrt(D_x*D_x + D_y*D_y);
        uL_D_x[i] = D_x/L_D_L[i];
        uL_D_y[i] = D_y/L_D_L[i];
        dL_D_L[i] = (dp_D_x[i1] - dp_D_x[i])*uL_D_x[i] + (dp_D_y[i1] - dp_D_y[i])*uL_D_y[i];
        L_V_L[i] = sqrt(V_x*V_x + V_y*V_y);
        uL_V_x[i] = V_x/L_V_L[i];
        uL_V_y[i] = V_y/L_V_L[i];
        dL_V_L[i] = (dp_V_x[i1] - dp_V_x[i])*uL_V_x[i] + (dp_V_y[i1] - dp_V_y[i])*uL_V_y[i];
    }
}


// Update the forces of the body given the current geometry and muscle activations

void WormBody::UpdateForces(int start, int end)
{
    int start2 = max(0,start-1), end2 = min(end+1,N_rods);
    double f_D_D_x[N_segments],f_D_D_y[N_segments],f_V_D_x[N_segments],f_V_D_y[N_segments];
    double f_D_L_x[N_segments],f_D_L_y[N_segments],f_V_L_x[N_segments],f_V_L_y[N_segments];
    double f_D_M_x[N_segments],f_D_M_y[N_segments],f_V_M_x[N_segments],f_V_M_y[N_segments];

    // Update individual forces
    for (int i = start2; i < end2 - 1; i++) {
        double temp;
        // Diagonal passive forces
        temp = kappa_D*(L_D0[i] - L_D_D[i]) - beta_D*dL_D_D[i];
        f_D_D_x[i] = temp*uD_D_x[i];
        f_D_D_y[i] = temp*uD_D_y[i];
        temp = kappa_D*(L_D0[i] - L_V_D[i]) - beta_D*dL_V_D[i];
        f_V_D_x[i] = temp*uD_V_x[i];
        f_V_D_y[i] = temp*uD_V_y[i];
        // Lateral passive forces
        temp = L_L0[i] - L_D_L[i];
        if (temp < 0) {double temp2 = temp*temp; temp += 16*temp2*temp2;}      // **** NOTE 2 ****
        temp = kappa_L*temp - beta_L*dL_D_L[i];
        f_D_L_x[i] = temp*uL_D_x[i];
        f_D_L_y[i] = temp*uL_D_y[i];
        temp =  L_L0[i] - L_V_L[i];
        if (temp < 0) {double temp2 = temp*temp; temp += 16*temp2*temp2;}      // **** NOTE 2 ****
        temp = kappa_L*temp - beta_L*dL_V_L[i];
        f_V_L_x[i] = temp*uL_V_x[i];
        f_V_L_y[i] = temp*uL_V_y[i];
        // Lateral active muscle forces
        temp = kappa_M0*A_D_M[i]*(L_L0[i] - A_D_M[i]*L_L0_minus_L_min[i] - L_D_L[i]) - beta_M0*A_D_M[i]*dL_D_L[i];
        f_D_M_x[i] = temp*uL_D_x[i];
        f_D_M_y[i] = temp*uL_D_y[i];
        temp = kappa_M0*A_V_M[i]*(L_L0[i] - A_V_M[i]*L_L0_minus_L_min[i] - L_V_L[i]) - beta_M0*A_V_M[i]*dL_V_L[i];
        f_V_M_x[i] = temp*uL_V_x[i];
        f_V_M_y[i] = temp*uL_V_y[i];
    }
    
    // Update total forces
    if (start2 == 0) {
        f_D_x[0] = -(f_D_D_x[0] + f_D_L_x[0] + f_D_M_x[0]);
        f_D_y[0] = -(f_D_D_y[0] + f_D_L_y[0] + f_D_M_y[0]);
        f_V_x[0] = -(f_V_D_x[0] + f_V_L_x[0] + f_V_M_x[0]);
        f_V_y[0] = -(f_V_D_y[0] + f_V_L_y[0] + f_V_M_y[0]);
    }
    for (int i = start2 + 1; i < end2 - 1; i++) {
        int i1 = i - 1;
        f_D_x[i] = f_V_D_x[i1] - f_D_D_x[i] + f_D_L_x[i1] - f_D_L_x[i] + f_D_M_x[i1] - f_D_M_x[i];
        f_D_y[i] = f_V_D_y[i1] - f_D_D_y[i] + f_D_L_y[i1] - f_D_L_y[i] + f_D_M_y[i1] - f_D_M_y[i];
        f_V_x[i] = f_D_D_x[i1] - f_V_D_x[i] + f_V_L_x[i1] - f_V_L_x[i] + f_V_M_x[i1] - f_V_M_x[i];
        f_V_y[i] = f_D_D_y[i1] - f_V_D_y[i] + f_V_L_y[i1] - f_V_L_y[i] + f_V_M_y[i1] - f_V_M_y[i];
    }
    if (end2 == N_rods) {
        f_D_x[N_segments] = f_V_D_x[N_segments-1] + f_D_L_x[N_segments-1] + f_D_M_x[N_segments-1];
        f_D_y[N_segments] = f_V_D_y[N_segments-1] + f_D_L_y[N_segments-1] + f_D_M_y[N_segments-1];
        f_V_x[N_segments] = f_D_D_x[N_segments-1] + f_V_L_x[N_segments-1] + f_V_M_x[N_segments-1];
        f_V_y[N_segments] = f_D_D_y[N_segments-1] + f_V_L_y[N_segments-1] + f_V_M_y[N_segments-1];
    }
}


// Update the residuals of the body states given the current forces and derivatives

void WormBody::UpdateResiduals(int start, int end)
{
    for (int i = start; i < end; i++) {
        int i3 = 3*i;
        double
        f_D_par = f_D_x[i]*sinPhi[i] - f_D_y[i]*cosPhi[i],
        f_D_perp = f_D_x[i]*cosPhi[i] + f_D_y[i]*sinPhi[i],
        f_V_par = f_V_x[i]*sinPhi[i] - f_V_y[i]*cosPhi[i],
        f_V_perp = f_V_x[i]*cosPhi[i] + f_V_y[i]*sinPhi[i],
        f_even_par = (f_V_par + f_D_par)/2,
        f_odd_par = (f_V_par - f_D_par)/2,
        V_CoM_perp = (f_D_perp + f_V_perp)/C_perp,
        V_CoM_par = 2*f_even_par/C_par;
            
        Residuals[i3] = dZ[i3] - (V_CoM_par*sinPhi[i] + V_CoM_perp*cosPhi[i]);
        Residuals[i3+1] = dZ[i3+1] - (-V_CoM_par*cosPhi[i] + V_CoM_perp*sinPhi[i]);
#ifdef BBC_STRICT
        Residuals[i3+2] = dZ[i3+2] - (f_odd_par/C_par)/(M_PI*2.0*R[i]);               // **** NOTE 3 ****
#else
        Residuals[i3+2] = dZ[i3+2] - 2*f_odd_par/(R[i]*C_par);                        // **** NOTE 3 ****
#endif
    }
}


// Numerically estimate the residual Jacobian matrix dF/dZ at the current state (Z, dZ) using finite differences
// Only computes the parts of dFdZ that can be nonzero
// This assumes that Residual already contains the residuals at (Z, dZ)
// Note: This perserves only Z, dZ and Residuals

void WormBody::NumericaldFdZ(double J[N_states][N_states])
{
    int j3, j33, lb, ub;
    double h, h1, temp, eps = sqrt(DBL_EPSILON);
    double Res[N_states];
    
    // Save the current residuals
    for (int i = 0; i < N_states; i++) Res[i] = Residuals[i];
    // Compute the finite difference approximiation
    for (int j = 0; j < N_states; j++) {
        // Perturb Z[j]
        temp = Z[j];
        h = eps*abs(temp);
        if (h == 0.0) h = eps;
        Z[j] = temp + h;
        // Trick to make sure h is exactly representable
        h = Z[j] - temp;
        // Only compute the the residuals that can change
        j3 = j/3; lb = max(0, j3-1); ub = min(j3+2, N_rods);
        F(lb, ub);
        // Restore Z[j]
        Z[j] = temp;
        // Update only the entries that can be nonzero
        h1 = 1/h;
        j33 = 3*j3;
        for (int i = max(0, j33-3); i < min(j33+6, N_states); i++) J[i][j] = h1*(Residuals[i] - Res[i]);
    }
    // Restore Residuals
    for (int i = 0; i < N_states; i++) Residuals[i] = Res[i];
}


// Numerically estimate the residual Jacobian matrix dF/dZ' at the current state (Z, dZ) using finite differences
// Only computes the parts of dFdZ' that can be nonzero
// This assumes that Residuals already contains the residuals at (Z, dZ)
// Note: This perserves only Z, dZ and Residuals

void WormBody::NumericaldFdZp(double J[N_states][N_states])
{
    int j3, j33, lb, ub;
    double h, h1, temp, eps = sqrt(DBL_EPSILON);
    double Res[N_states];
    
    // Save the current residuals
    for (int i = 0; i < N_states; i++) Res[i] = Residuals[i];
    // Compute the finite difference approximiation
    for (int j = 0; j < N_states; j++) {
        // Perturb dZ[j]
        temp = dZ[j];
        h = eps*abs(temp);
        if (h == 0.0) h = eps;
        dZ[j] = temp + h;
        // Trick to make sure h is exactly representable
        h = dZ[j] - temp;
        // Only compute the the residuals that can change
        j3 = j/3; lb = max(0, j3-1); ub = min(j3+2, N_rods);
        F(lb,ub);
        // Restore dZ[j]
        dZ[j] = temp;
        // Update only the entries that can be nonzero
        h1 = 1/h;
        j33 = 3*j3;
        for (int i = max(0, j33-3); i < min(j33+6, N_states); i++) J[i][j] = h1*(Residuals[i] - Res[i]);
    }
    // Restore Residuals
    for (int i = 0; i < N_states; i++) Residuals[i] = Res[i];
}


// Solve the linear system M . Z = B using Gaussian elimination w/o pivoting
// This is faster than LU decomposition, but may be less stable
// This code takes advantage of the banded structure of M
// Modifies M and Z, but not B

void WormBody::LinearSolve(double M[N_states][N_states], double B[], double Z[])
{
    const int bw2 = 5;
    int i, j, k;
    
    // Copy B to Z
    for (i = 0; i < N_states; i++) Z[i] = B[i];
    // Peform Gaussian elimination (w/o pivoting) and forward substitution w/in the band
    for (k = 0; k < N_states - 1; k++)
        for (i = k + 1; i < min(k+bw2+1, N_states); i++) {
            double temp = M[i][k] /= M[k][k];
            for (j = k + 1; j < min(k+bw2+1, N_states); j++)
                M[i][j] -= temp*M[k][j];
            Z[i] -= temp*Z[k];
        }
    // Perform backsubstitution w/in the band
    for (i = N_states - 1; i >= 0; i--) {
        double sum = Z[i];
        for (int j = i + 1; j < min(i+6, N_states); j++) sum -= M[i][j]*Z[j];
        Z[i] = sum/M[i][i];
    }
}

