#include <fstream>
#include <vector>

#include "Collide.h"

std::vector<CollisionObject> load_objects()
{
	std::vector<CollisionObject> CollObjs = std::vector<CollisionObject>();

    // open file
    std::ifstream objfile(COLLIDE_FILE);
    if (!objfile.is_open())
    {
        exit(EXIT_FAILURE);
    }

    // initialize temp variables
    double bound_min_x;
    double bound_min_y;

    double bound_max_x;
    double bound_max_y;

    double fvec_x;
    double fvec_y;

    // loop
    while(
        objfile 
        >> bound_min_x >> bound_min_y 
        >> bound_max_x >> bound_max_y 
        >> fvec_x >> fvec_y
    ){
        CollisionObject tempObj;
        tempObj.bound_min_x = bound_min_x; 
        tempObj.bound_min_y = bound_min_y;
        tempObj.bound_max_x = bound_max_x; 
        tempObj.bound_max_y = bound_max_y; 
        tempObj.fvec_x = fvec_x; 
        tempObj.fvec_y = fvec_y;
        
        // store data
        CollObjs.push_back(tempObj);
    }

    // close file
    objfile.close();

	return CollObjs;
}