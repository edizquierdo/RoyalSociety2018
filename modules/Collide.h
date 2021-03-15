#ifndef COLLIDE_H
#define COLLIDE_H

#include <vector>
#include <fstream>
// #include <stdlib.h> 
// #include <cmath>

#define COLLIDE_FILE "data/collision_objs.tsv"


struct CollisionObject
{
	double bound_min_x;
	double bound_min_y;

	double bound_max_x;
	double bound_max_y;

	double fvec_x;
	double fvec_y;
};

std::vector<CollisionObject> load_objects();

#endif

/* 
collision code originally from:
@article{Boyle_Berri_Cohen_2012, 
 	title={Gait Modulation in C. elegans: An Integrated Neuromechanical Model}, 
	volume={6}, 
	ISSN={1662-5188}, 
	url={https://www.frontiersin.org/articles/10.3389/fncom.2012.00010/full#h8}, 
 	DOI={10.3389/fncom.2012.00010}, 
	journal={Frontiers in Computational Neuroscience}, 
	publisher={Frontiers}, 
	author={Boyle, Jordan Hylke and Berri, Stefano and Cohen, Netta}, 
	year={2012}
modified by github.com/mivanit
} */


/*
	// If using objects, check for object collisions and calculate associated forces
	if(N_objects > 0){
		// REVIEW: unclear what this does
  		realtype P_x,P_y,Distance,magF,D_scale,magF_P1,magF_P2;
		// reset contact force
  		ContactForce = 0;

		// for each side of the solid rods that make up the worm body
  		for(int i = 0; i < NBAR; ++i){
			for(int j = 0; j < 2; ++j){
				// First ensure they contain zeros
				F_object[i][j][0] = 0;
				F_object[i][j][1] = 0;
				P_x = term[i][j][0];
				P_y = term[i][j][1];

				// Now check proximity to each object
				for(int k = 0; k < N_objects; ++k){
					// simple bounding box check
					if(
							(P_x < (Objects[k][0]+Objects[k][2]))
							&& (P_x > (Objects[k][0]-Objects[k][2]))
							&& (P_y < (Objects[k][1]+Objects[k][2]))
							&& (P_y > (Objects[k][1]-Objects[k][2]))
						){

						// distance to the center of the object. this will need to be changed for non-disc objects
						dx = P_x - Objects[k][0];
						dy = P_y - Objects[k][1];
						Distance = sqrt(pow(dx,2) + pow(dy,2));
						D_scale = 0.01*Objects[k][2];

						// TODO: instead of checking distance, check inside collision box (see existing python code)
						if(Distance < Objects[k][2]){
							magF = (
								k_Object * (Objects[k][2] - Distance) 
								+ D_scale*k_Object * (pow((Objects[k][2] - Distance) / D_scale,2))
							);
							// this updates the object forces based on how deep into the object we are
							// `i` is the segment, and `j` is the ventral/dorsal side
							F_object[i][j][0] += (dx/Distance)*magF;
							F_object[i][j][1] += (dy/Distance)*magF;
							ContactForce += magF;
						}
					}			
				}
			}
  		}
	}

  	// Add up force contributions for each D/V point
  	F_term[0][0][0] = -F_H[0][0]*Dir[0][0][0] - F_D[0][0]*Dir_D[0][0][0] + F_object[0][0][0];
  	F_term[0][0][1] = -F_H[0][0]*Dir[0][0][1] - F_D[0][0]*Dir_D[0][0][1] + F_object[0][0][1];

  	F_term[0][1][0] = -F_H[0][1]*Dir[0][1][0] - F_D[0][1]*Dir_D[0][1][0] + F_object[0][1][0];
  	F_term[0][1][1] = -F_H[0][1]*Dir[0][1][1] - F_D[0][1]*Dir_D[0][1][1] + F_object[0][1][1];

	// loop over each segment adding force
	// TODO: this needs to call the data structure used in the Izquierdo code
  	for(int i = 1; i < NSEG; ++i){
		int i_minus_1 = i-1;

		F_term[i][0][0] = F_H[i_minus_1][0]*Dir[i_minus_1][0][0] - F_H[i][0]*Dir[i][0][0] + F_D[i_minus_1][1]*Dir_D[i_minus_1][1][0] - F_D[i][0]*Dir_D[i][0][0] + F_object[i][0][0];
		F_term[i][0][1] = F_H[i_minus_1][0]*Dir[i_minus_1][0][1] - F_H[i][0]*Dir[i][0][1] + F_D[i_minus_1][1]*Dir_D[i_minus_1][1][1] - F_D[i][0]*Dir_D[i][0][1] + F_object[i][0][1];

		F_term[i][1][0] = F_H[i_minus_1][1]*Dir[i_minus_1][1][0] - F_H[i][1]*Dir[i][1][0] + F_D[i_minus_1][0]*Dir_D[i_minus_1][0][0] - F_D[i][1]*Dir_D[i][1][0] + F_object[i][1][0];
		F_term[i][1][1] = F_H[i_minus_1][1]*Dir[i_minus_1][1][1] - F_H[i][1]*Dir[i][1][1] + F_D[i_minus_1][0]*Dir_D[i_minus_1][0][1] - F_D[i][1]*Dir_D[i][1][1] + F_object[i][1][1];
  	}

  	F_term[NSEG][0][0] = F_H[NSEG_MINUS_1][0]*Dir[NSEG_MINUS_1][0][0] + F_D[NSEG_MINUS_1][1]*Dir_D[NSEG_MINUS_1][1][0] + F_object[NSEG][0][0];
  	F_term[NSEG][0][1] = F_H[NSEG_MINUS_1][0]*Dir[NSEG_MINUS_1][0][1] + F_D[NSEG_MINUS_1][1]*Dir_D[NSEG_MINUS_1][1][1] + F_object[NSEG][0][1];

  	F_term[NSEG][1][0] = F_H[NSEG_MINUS_1][1]*Dir[NSEG_MINUS_1][1][0] + F_D[NSEG_MINUS_1][0]*Dir_D[NSEG_MINUS_1][0][0] + F_object[NSEG][1][0];
  	F_term[NSEG][1][1] = F_H[NSEG_MINUS_1][1]*Dir[NSEG_MINUS_1][1][1] + F_D[NSEG_MINUS_1][0]*Dir_D[NSEG_MINUS_1][0][1] + F_object[NSEG][1][1];
  
  	// Convert net forces on D/V points to force and torque	acting on rod CoM
  	for(int i = 0; i < NBAR; ++i){
		realtype cos_thi = cos(CoM[i][2]);
		realtype sin_thi = sin(CoM[i][2]);
		for(int j = 0; j < 2; ++j){			
			F_term_rotated[i][j][0] = F_term[i][j][0]*cos_thi + F_term[i][j][1]*sin_thi;	// This is Fperp
			F_term_rotated[i][j][1] = F_term[i][j][0]*sin_thi - F_term[i][j][1]*cos_thi;    // THis is Fparallel
		}

		V_CoM_rotated[i][0] = (F_term_rotated[i][0][0] + F_term_rotated[i][1][0])/CN[i];

		F_even = (F_term_rotated[i][0][1] + F_term_rotated[i][1][1]);	//Took out the /2
		F_odd = (F_term_rotated[i][1][1] - F_term_rotated[i][0][1])/RCONST(2.0);	

		V_CoM_rotated[i][1] = (F_even)/CL[i];				//Allowing me to take out *2
		V_CoM[i][2] = (F_odd/CL[i])/(M_PI*2.0*R[i]);

		V_CoM[i][0] = V_CoM_rotated[i][0]*cos_thi + V_CoM_rotated[i][1]*sin_thi;
		V_CoM[i][1] = V_CoM_rotated[i][0]*sin_thi - V_CoM_rotated[i][1]*cos_thi;
	
		int three_i = i*3;

		rval[three_i] = V_CoM[i][0] - ypval[three_i];
		rval[three_i+1] = V_CoM[i][1] - ypval[three_i+1];
		rval[three_i+2] = V_CoM[i][2] - ypval[three_i+2];
  	}
*/

