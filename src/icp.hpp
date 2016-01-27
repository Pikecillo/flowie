#ifndef __ICP_HPP__
#define __ICP_HPP__

#include <vector>

#include "matrix_types.hpp"
#include "vector_types.hpp"

Mat3x4d affine_registration(std::vector<Vector3d> source,
			    std::vector<Vector3d> target);
Mat3x4d icp(std::vector<Vector3d> source,
	    std::vector<Vector3d> target,
	    int max_iter);

#endif
