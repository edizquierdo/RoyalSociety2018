#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Collide.h"

// #define COLLIDE_DEBUG

#ifdef COLLIDE_DEBUG
	#include <iostream>
#endif

VecXY get_displacement(VecXY a, VecXY b)
{
	return VecXY(
		a.x - b.x,
		a.y - b.y
	);
}

double dist(VecXY a, VecXY b)
{
	return pow((
		pow(a.x - b.x, 2.0)
		+ pow(a.y - b.y, 2.0)
	), 0.5);
}

std::vector<CollisionObject> load_objects(std::string collide_file)
{
	std::vector<CollisionObject> CollObjs = std::vector<CollisionObject>();

    // open file
    std::ifstream objfile(collide_file);
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

			// store data
			CollObjs.push_back(tempObj);
		}
    }

    // close file
    objfile.close();

	#ifdef COLLIDE_DEBUG
		for (CollisionObject obj : CollObjs)
		{
			std::cout << obj.coll_type << "," << obj.bound_max_x << "," << obj.bound_min_x << "," << obj.force << std::endl;
		}
	#endif

	return CollObjs;
}


void save_objects(std::string collide_file, std::vector<CollisionObject> & CollObjs)
{
    // open file
    std::ofstream objfile(collide_file);
    if (!objfile.is_open() || !objfile.good())
    {
        exit(EXIT_FAILURE);
    }

	PRINTF_DEBUG("    >> elements in CollObjs vec: %d\n", CollObjs.size())

	for (CollisionObject obj : CollObjs)
	{
		if (obj.coll_type == Box_Ax)
		{
			objfile 
				<< "Box_Ax"
				<< "\t" << obj.bound_min_x << "\t" << obj.bound_min_y 
				<< "\t" << obj.bound_max_x << "\t" << obj.bound_max_y 
				<< "\t" << obj.fvec_x << "\t" << obj.fvec_y
				<< std::endl;
		}
		else if (obj.coll_type == Disc)
		{
			objfile 
				<< "Disc"
				<< "\t" << obj.bound_min_x << "\t" << obj.bound_min_y 
				<< "\t" << obj.bound_max_x << "\t" << obj.bound_max_y 
				<< "\t" << obj.centerpos_x << "\t" << obj.centerpos_y
				<< "\t" << obj.force
				<< "\t" << obj.radius_inner << "\t" << obj.radius_outer
				<< "\t" << obj.angle_min << "\t" << obj.angle_max
				<< std::endl;
		}
		else
		{
			objfile << "NULL\n";
		}
    }

    // close file
	objfile.flush();
    objfile.close();
}


VecXY do_collide(CollisionObject obj, VecXY pos)
{
	// forces on elements
	if (obj.coll_type == Box_Ax)
	{			
		if (
			(pos.x > obj.bound_min_x)
			&& (pos.x < obj.bound_max_x)
			&& (pos.y > obj.bound_min_y)
			&& (pos.y < obj.bound_max_y)
		){
			return VecXY(obj.fvec_x, obj.fvec_y);
		}
	}
	else if (obj.coll_type == Disc)
	{
		VecXY disc_center = VecXY(obj.centerpos_x, obj.centerpos_y);
		VecXY offset = get_displacement(pos, disc_center);
		double offset_mag = offset.mag();
		// compare to radius
		if (
			( offset_mag > obj.radius_inner) 
			&& (offset_mag < obj.radius_outer)
		)
		{
			// check wedge angles
			double angle = atan2(offset.y, offset.x);			
			if ( !( (obj.angle_min > angle) && (angle > obj.angle_max) ) )
			{
				// get the collision vector by normalizing and scaling offset
				offset.scale(obj.force / offset_mag);
				return offset;
			}
		}
	}

	// if no collisions found, return zero vec
	return VecXY();
}
