/*
 * University of Houston
 * Mario Rincon-Nigro. May 2013.
 */

#ifndef __GRAPHICS_HPP__
#define __GRAPHICS_HPP__

#include <GL/gl.h>

#include "color_scheme.hpp"
#include "geometry.hpp"
#include "histogram.hpp"
#include "vector_types.hpp"

#include "mesh/triangle_soup.hpp"

template <class S_>
void draw_box(const AABB<S_, 3> &aabb) {
  Vector<S_, 3> corners[] = { aabb.min, aabb.max };
  int order[] = { 0, 1, 3, 2, 4, 5, 7, 6 };

  // Draw rectangles for two of the opposite sides
  glBegin(GL_LINE_LOOP);
  for(int i = 0; i < 8; i++) {
    int x = (order[i] & 1), y = (order[i] & 2) >> 1, z = (order[i] & 4) >> 2;
    glVertex3f(corners[x].x(), corners[y].y(), corners[z].z());

    if(i == 3) {
      glEnd();
      glBegin(GL_LINE_LOOP);
    }
  }
  glEnd();

  // Connect the rectangles
  glBegin(GL_LINES);
  for(int i = 0; i < 4; i++) {
    int x = (i & 1), y = (i & 2) >> 1;
    glVertex3f(corners[x].x(), corners[y].y(), corners[0].z());
    glVertex3f(corners[x].x(), corners[y].y(), corners[1].z());
  }
  glEnd();
}

template <class S_>
void draw_wireframe(PolygonSoupOBJ<S_> &obj_mesh) {
  // For each polygon
  for(unsigned int i = 0; i < obj_mesh.polygons.size(); i++) {
    FaceOBJ tri = obj_mesh.polygons[i];

    glBegin(GL_LINE_LOOP);
    for(unsigned int i = 0; i < tri.size(); i++) {
      Vector<S_, 3> v = obj_mesh.positions[tri[i][0]];
      glVertex3f(v.x(), v.y(), v.z());
    }
    glEnd();
  }
}
 
template <class S_>
void draw_points(const std::vector<Vector<S_, 2> > &points,
		 const Color3f &color) {
  glPointSize(2.0);
  glColor3f(color[0], color[1], color[2]);

  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < points.size(); i++) {
    Vector<S_, 2> p = points[i];
    glVertex2f(p.x(), p.y());
  }
  glEnd();
}

template <class S_>
void draw_points(const std::vector<Vector<S_, 3> > &points,
                 const Color3f &color) {
  glPointSize(2.0);
  glColor3f(color[0], color[1], color[2]);

  glBegin(GL_POINTS);
  for(unsigned int i = 0; i < points.size(); i++) {
    Vector<S_, 3> p = points[i];
    glVertex3f(p.x(), p.y(), p.z());
  }
  glEnd();
}

template <class S_>
void draw_vector(const Vector<S_, 2> &origin,
                 const Vector<S_, 2> &direction,
                 bool normalize) {
  S_ scale_b = 0.2, scale_h = 0.25;
  Vector<S_, 2> dir = direction;
  
  if(normalize)
    dir = dir.normalize() * 0.3;

  // Draw the arrow line
  glPushMatrix();
  glTranslatef(origin[0], origin[1], 0);
  glScalef(scale_b, scale_b, 1);
  glBegin(GL_LINES);
  glVertex2f(0.0f, 0.0f);
  glVertex2f(dir[0], dir[1]);
  glEnd();
  glPopMatrix();
  
  // Draw arrow head
  glPushMatrix();
  glTranslatef(origin[0], origin[1], 0);
  glScalef(scale_b, scale_b, 1);
  glTranslatef(dir[0], dir[1], 0);
  glRotatef(atan2(dir[1], dir[0]) * 360 / (2 * M_PI), 0, 0, 1);
  glScalef(scale_h, scale_h, 1);
  glBegin(GL_TRIANGLES);
  glVertex2f(0, 0);
  glVertex2f(-0.35, 0.12);
  glVertex2f(-0.35, -0.12);
  glEnd();
  glPopMatrix();
}

template <class S_>
void draw_vector_field(const VectorField<S_, 2> &field,
                       bool normalize) {
  Histogram<S_> histogram(512);
  ScalarField<S_, 2> result(field.getConf());
  
  for(int i = 0; i < result.size(); i++) {
    result.at(i) = field.at(i).length();
  }

  build_histogram(result, histogram);
  
  for(int i = 0; i < field.size(); i++) {
    Vector2d v = field.at(i);
    Color3f rgb = ColorScheme::getColor(v.length(), histogram,
    					ColorScheme::BWR, true);
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    draw_vector(field.coordinateAt(i), field.at(i), normalize);
  }
}

template <class S_, int dim_>
void draw_isocontour
(const std::vector<std::pair<Vector<S_, dim_>,
 Vector<S_, dim_> > > &isocontour, const Color3f &color) {
  // Draw line segments
  glColor3f(color[0], color[1], color[2]);
  
  glBegin(GL_LINES);
  for(unsigned int i = 0; i < isocontour.size(); i++) {
    Vector<S_, 2> v0 = isocontour[i].first, v1 = isocontour[i].second;
    glVertex3f(v0.x(), v0.y(), 0.0f);
    glVertex3f(v1.x(), v1.y(), 0.0f);
  }
  glEnd();
}

template <class S_>
 void draw_scalar_field_grid(const ScalarField<S_, 2> &field,
                             S_ threshold,
                             bool signed_regions) {
  Histogram<S_> histogram(512);
  build_histogram(field, histogram);

  glBegin(GL_POINTS);
  
  for(int i = 0; i < field.size(); i++) {
    Vector<S_, 2> point = field.coordinateAt(i);
    Color3f rgb = ColorScheme::getColor(field.at(i), histogram,
                                        ColorScheme::BWR, true);
    
    if(signed_regions) {
      if(field.at(i) >= 0)
        rgb = Color3f(0.0, 0.0, 1.0);
      else
        rgb = Color3f(1.0, 0.0, 0.0);
    }

    if(field.at(i) <= threshold) {
      glColor3f(rgb[0], rgb[1], rgb[2]);
      glVertex3d(point.x(), point.y(), 0.0f);
    }
  }
  
  glEnd();
}

template <class S_>
 void draw_scalar_field_grid(const ScalarField<S_, 3> &field,
                             S_ threshold,
                             bool signed_regions) {
  Histogram<S_> histogram(512);
  build_histogram(field, histogram);

  glBegin(GL_POINTS);
  
  for(int i = 0; i < field.size(); i++) {
    Vector<S_, 3> point = field.coordinateAt(i);
    Color3f rgb = ColorScheme::getColor(field.at(i), histogram,
                                        ColorScheme::BWR, true);
    
    if(signed_regions) {
      if(field.at(i) >= 0)
        rgb = Color3f(0.0, 0.0, 1.0);
      else
        rgb = Color3f(1.0, 0.0, 0.0);
    }

    if(field.at(i) <= threshold) {
      glColor3f(rgb[0], rgb[1], rgb[2]);
      glVertex3d(point.x(), point.y(), point.z());
    }
  }
  
  glEnd();
}

template <class S_>
void draw_scalar_field_solid(const ScalarField<S_, 2> &slice,
			     int orientation, S_ free_coord,
			     Histogram<S_> &histogram,
			     bool signed_regions) {
  int y_dim = slice.getConf().numCells.y();
  int x_dim = slice.getConf().numCells.x();

  std::cout << "-----" << std::endl;
  for(int i = 0; i < y_dim; i++) {
    for(int j = 0; j < x_dim; j++) {
      std::cout << i << " " << j << " " << slice.at(i, j) << std::endl;

      if(i >= y_dim - 1) continue;
      if(j >= x_dim - 1) continue;

      Vector<S_, 2> p[] = {
	slice.coordinateAt(i, j),
	slice.coordinateAt(i, j + 1),
	slice.coordinateAt(i + 1, j),
	slice.coordinateAt(i + 1, j + 1)
      };
      Color3f rgb[] = {
	ColorScheme::getColor(slice.at(i, j), histogram,
			      ColorScheme::BWR, true),
	ColorScheme::getColor(slice.at(i, j + 1), histogram,
			      ColorScheme::BWR, true),
	ColorScheme::getColor(slice.at(i + 1, j), histogram,
			      ColorScheme::BWR, true),
	ColorScheme::getColor(slice.at(i + 1, j + 1), histogram,
			      ColorScheme::BWR, true)
      };

      glBegin(GL_TRIANGLE_STRIP);
      for(int i = 0; i < 4; i++) {
	if(signed_regions) {
	  if(slice.at(i, j) >= 0)
	    rgb[i] = Color3f(0.0, 0.0, 1.0);
	  else
	    rgb[i] = Color3f(1.0, 0.0, 0.0);
	}
	
	glColor3f(rgb[i][0], rgb[i][1], rgb[i][2]);
	if(orientation == XY)
	  glVertex3d(p[i].x(), p[i].y(), free_coord);
	else if(orientation == XZ)
	  glVertex3d(p[i].x(), free_coord, p[i].y());
	else
	  glVertex3d(free_coord, p[i].x(), p[i].y());
      }
      glEnd();
    }
  }
}

#endif
