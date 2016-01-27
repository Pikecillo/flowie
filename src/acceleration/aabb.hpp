/*
 * University of Houston
 * Mario Rincon-Nigro. July 2013.
 */

#ifndef __AABB_HPP__
#define __AABB_HPP__

#include <assert.h>

#include "geometry.hpp"
#include "vector_types.hpp"

template <class T_, int dim_>
struct AABB {
  Vector<T_, dim_> max;
  Vector<T_, dim_> min;
  
  AABB() {
    T_ T_max = std::numeric_limits<T_>::max();
    max = Vector<T_, dim_>(-T_max);
    min = Vector<T_, dim_>(T_max);
  }
  
  AABB(const AABB<T_, dim_> &other) {
    (*this) = other;
  }
  
  AABB(const Vector<T_, dim_> &omin, const Vector<T_, dim_> &omax) {
    assert(omin <= omax);
    min = omin;
    max = omax;
  }
  
  const AABB<T_, dim_> &operator=(const AABB<T_, dim_> &other) {
    min = other.min;
    max = other.max;
    
    return (*this);
  }
  
  void grow(const Vector<T_, dim_> &point) {
    min = min.min(point);
    max = max.max(point);
  }
  
  void grow(const Vector<T_, dim_> &p1,
	    const Vector<T_, dim_> &p2,
	    const Vector<T_, dim_> &p3) {
    grow(p1); grow(p2); grow(p3);
  }
  
  void grow(const std::vector<Vector<T_, dim_> > &points) {
    for(unsigned int i = 0; i < points.size(); i++)
      grow(points[i]);
  }

  void join(const AABB<T_, dim_> &other) {
    grow(other.min);
    grow(other.max);
  }
  
  void intersect(const AABB<T_, dim_> &other) {
    min = min.max(other.min);
    max = max.min(other.max);
  }
  
  T_ area() const {
    Vector<T_, dim_> diagonal = max - min;
    T_ area = (T_)0;
    
    for(int i = 0; i < diagonal.size(); i++)
      for(int j = i + 1; j < diagonal.size(); j++)
	area += (diagonal[i] * diagonal[j]);

    return 2 * area;
  }
  
  T_ volume() const {
    return (max - min).volume();
  }
  
  Vector<T_, 3> centroid() const {
    return (T_)0.5 * (max + min);
  }
};

/*
 * Minimum distance from point to axis-aligned bounding box
 */
template <class T_>
T_ min_distance(const AABB<T_, 3> &aabb, const Vector<T_, 3> &p) {  
  enum { MIN, MAX, MID };
  Vector<T_, 3> max = aabb.max, min = aabb.min;
  Vector<T_, 3> corners[8] = {
    Vector<T_, 3>(min.x(), min.y(), min.y()),
    Vector<T_, 3>(max.x(), min.y(), min.y()),
    Vector<T_, 3>(min.x(), max.y(), min.y()),
    Vector<T_, 3>(max.x(), max.y(), min.y()),
    Vector<T_, 3>(min.x(), min.y(), max.y()),
    Vector<T_, 3>(max.x(), min.y(), max.y()),
    Vector<T_, 3>(min.x(), max.y(), max.y()),
    Vector<T_, 3>(max.x(), max.y(), max.y())
    };
  Vector3i region;

  // Identify region
  region[0] = (p.x() < min.x() ? MIN : (p.x() > max.x() ? MAX : MID));
  region[1] = (p.y() < min.y() ? MIN : (p.y() > max.y() ? MAX : MID));
  region[2] = (p.z() < min.z() ? MIN : (p.z() > max.z() ? MAX : MID));

  // Check for corners
  if(region == Vector3i(MIN, MIN, MIN)) return (corners[0] - p).length();
  if(region == Vector3i(MAX, MIN, MIN)) return (corners[1] - p).length();
  if(region == Vector3i(MIN, MAX, MIN)) return (corners[2] - p).length();
  if(region == Vector3i(MAX, MAX, MIN)) return (corners[3] - p).length();
  if(region == Vector3i(MIN, MIN, MAX)) return (corners[4] - p).length();
  if(region == Vector3i(MAX, MIN, MAX)) return (corners[5] - p).length();
  if(region == Vector3i(MIN, MAX, MAX)) return (corners[6] - p).length();
  if(region == Vector3i(MAX, MAX, MAX)) return (corners[7] - p).length();
  
  // Check for sides for points outside
  T_ dsides[] = { min.z() - p.z(), p.z() - max.z(),
		  min.y() - p.y(), p.y() - max.y(), 
		  min.x() - p.x(), p.x() - max.x() };
  if(region == Vector3i(MID, MID, MIN)) return dsides[0];
  if(region == Vector3i(MID, MID, MAX)) return dsides[1];
  if(region == Vector3i(MID, MIN, MID)) return dsides[2];
  if(region == Vector3i(MID, MAX, MID)) return dsides[3];
  if(region == Vector3i(MIN, MID, MID)) return dsides[4];
  if(region == Vector3i(MAX, MID, MID)) return dsides[5];

  // Check sides for points inside
  if(region == Vector3i(MID, MID, MID)) {
    T_ dmin = -dsides[0];
      for(int i = 1; i < 6; i++)
	if(dmin > -dsides[i]) dmin = -dsides[i];

    return dmin;
  }

  // Check for edges
  if(region == Vector3i(MID, MAX, MIN))
    return min_distance(Segment<T_>(corners[2], corners[3]), p);
  if(region == Vector3i(MID, MAX, MAX))
    return min_distance(Segment<T_>(corners[6], corners[7]), p);
  if(region == Vector3i(MID, MIN, MAX))
    return min_distance(Segment<T_>(corners[4], corners[5]), p);
  if(region == Vector3i(MID, MIN, MIN))
    return min_distance(Segment<T_>(corners[0], corners[1]), p);
  if(region == Vector3i(MIN, MID, MIN))
    return min_distance(Segment<T_>(corners[0], corners[2]), p);
  if(region == Vector3i(MAX, MID, MIN))
    return min_distance(Segment<T_>(corners[1], corners[3]), p);
  if(region == Vector3i(MAX, MID, MAX))
    return min_distance(Segment<T_>(corners[5], corners[7]), p);
  if(region == Vector3i(MIN, MID, MAX))
    return min_distance(Segment<T_>(corners[4], corners[6]), p);
  if(region == Vector3i(MIN, MAX, MID))
    return min_distance(Segment<T_>(corners[2], corners[6]), p);
  if(region == Vector3i(MAX, MAX, MID))
    return min_distance(Segment<T_>(corners[3], corners[7]), p);
  if(region == Vector3i(MAX, MIN, MID))
    return min_distance(Segment<T_>(corners[1], corners[5]), p);
  if(region == Vector3i(MIN, MIN, MID))
    return min_distance(Segment<T_>(corners[0], corners[4]), p);

  assert(false);
}

#endif
