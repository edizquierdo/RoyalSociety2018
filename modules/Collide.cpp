#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Collide.h"

std::vector<CollisionObject> load_objects()
{
	std::vector<CollisionObject> CollObjs = std::vector<CollisionObject>();

    // open file
    std::ifstream objfile(COLLIDE_FILE);
    if (!objfile.is_open() || !objfile.good())
    {
        exit(EXIT_FAILURE);
    }

    // initialize temp variables
	
	// stores type
	std::string raw_line;
	std::string str_colltype;

    // loop
    while(getline(objfile, raw_line))
	{
		std::istringstream liness(raw_line);
		liness >> str_colltype;

		if (str_colltype == "Box_Ax")
		{
			CollisionObject tempObj;
			tempObj.coll_type = Box_Ax;

			liness 
				>> tempObj.bound_min_x >> tempObj.bound_min_y 
				>> tempObj.bound_max_x >> tempObj.bound_max_y 
				>> tempObj.fvec_x >> tempObj.fvec_y;

			// store data
			CollObjs.push_back(tempObj);
		}
		else if (str_colltype == "Disc")
		{
			CollisionObject tempObj;
			tempObj.coll_type = Disc;

			liness 
				>> tempObj.bound_min_x >> tempObj.bound_min_y 
				>> tempObj.bound_max_x >> tempObj.bound_max_y 
				>> tempObj.centerpos_x >> tempObj.centerpos_y
				>> tempObj.force
				>> tempObj.radius_inner >> tempObj.radius_outer
				>> tempObj.angle_min >> tempObj.angle_max;
		}
    }

    // close file
    objfile.close();

	return CollObjs;
}