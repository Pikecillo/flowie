/*
 * color_scheme.hpp
 *
 *  Created on: Mar 27, 2013
 *      Author: mario
 */

#ifndef __COLOR_SCHEME_HPP__
#define __COLOR_SCHEME_HPP__

#include <math.h>

#include <algorithm>

#include "histogram.hpp"
#include "vector_types.hpp"

typedef Vector3f Color3f;

class ColorScheme {
public:
  enum Scheme { RAINBOW, BWR, HEAT, GRAY };

public:
  /*
   * Routine to convert HSV to RGB
   *
   * Reference:  Foley, van Dam, Feiner, Hughes,
   * "Computer Graphics Principles and Practices,"
   * Additon-Wesley, 1990, pp 592-593.
   */
  static Color3f hsvToRgb(const Color3f &hsv) {
    double h, s, v;	  // hue, sat, value
    double r, g, b;	  // red, green, blue
    double i, f, p, q, t;  // interim values
    Color3f rgb;

    // Guarantee valid input:
    h = hsv[0] / 60.0;
    
    while(h >= 6.0) h -= 6.;
    while(h <  0.0) h += 6.;
    
    s = hsv[1];
    if(s < 0.0) s = 0.0;
    if(s > 1.0) s = 1.0;
    
    v = hsv[2];
    if(v < 0.0) v = 0.0;
    if(v > 1.0) v = 1.0;
    
    // If sat is 0, then is a gray:
    if(s == 0.0) {
      rgb[0] = rgb[1] = rgb[2] = v;
      return rgb;
    }
    
    // Get an rgb from the hue itself:
    i = floor(h);
    f = h - i;
    p = v * (1.0 - s);
    q = v * (1.0 - s * f );
    t = v * (1.0 - (s * (1.0 - f)));
    
    switch((int) i) {
    case 0:
      r = v; g = t; b = p;
      break;
    case 1:
      r = q; g = v; b = p;
      break;
    case 2:
      r = p; g = v; b = t;
      break;
    case 3:
      r = p; g = q; b = v;
      break;
    case 4:
      r = t; g = p; b = v;
      break;
    case 5:
      r = v; g = p; b = q;
      break;
    }
    
    rgb[0] = r; rgb[1] = g; rgb[2] = b;

    return rgb;
  }
  
  /*
   * Reference: http://www.cs.rit.edu/~ncs/color/t_convert.html
   */
  static Color3f rgbToHsv(const Color3f &rgb) {
    double min, max, delta;
    Color3f hsv;

    min = std::min(rgb[0], std::min(rgb[1], rgb[2]));
    max = std::max(rgb[0], std::max(rgb[1], rgb[2]));
    hsv[2] = max;
    
    delta = max - min;
    
    if(max != 0)
      hsv[1] = delta / max;
    else {
      // If r = g = b = 0 then s = 0, v is undefined
      hsv[1] = 0;
      hsv[0] = -1;
      return rgb;
    }
    
    if(rgb[0] == max)
      hsv[0] = (rgb[1] - rgb[2]) / delta; // Between yellow & magenta
    else if(rgb[2] == max)
      hsv[0] = 2 + (rgb[2] - rgb[0]) / delta; // Between cyan & yellow
    else
      hsv[0] = 4 + (rgb[0] - rgb[1]) / delta; // Between magenta & cyan
    
    hsv[0] *= 60; // degrees

    if(hsv[0] < 0)
      hsv[0] += 360;

    return rgb;
  }

  // Transfer functions
  static Color3f rainbowMap(double s) {
    return hsvToRgb(Color3f(240.0f * (1.0f - s), 1.0, 1.0));
  }
  
  static Color3f heatMap(double s) {
    double one_third = 1.0 / 3.0;
    Color3f rgb;

    rgb[0] = rgb[1] = rgb[2] = 0.0;
    
    if(s < one_third) {
      rgb[0] = 3.0 * s;
    } else if(s < 2.0 * one_third) {
      rgb[0] = 1.0f;
      rgb[1] = 3.0 * s - 1.0f;
    } else {
      rgb[0] = rgb[1] = 1.0f;
      rgb[2] = 3.0 * s - 2.0f;
    }

    return rgb;
  }
  
  static Color3f bwrMap(double s) {
    Color3f rgb;

    if(s < 0.5) {
      rgb[0] = rgb[1] = 2.0 * s;
      rgb[2] = 1.0;
    } else {
      rgb[0] = 1.0;
      rgb[1] = rgb[2] = 2.0 * (1.0 - s);
    }

    return rgb;
  }
  
  static Color3f grayMap(double s) {
    return Color3f(s, s, s);
  }
  
  static Color3f getColor(double value, const Histogram<float> &histogram,
                          Scheme scheme=RAINBOW, bool adjust_range=false,
                          bool adjust_hsv=false) {
    return getColor<float>(value, histogram, scheme,
                           adjust_range, adjust_hsv);
  }

  static Color3f getColor(double value, const Histogram<double> &histogram,
                          Scheme scheme=RAINBOW, bool adjust_range=false,
                          bool adjust_hsv=false) {
    return getColor<double>(value, histogram, scheme,
                            adjust_range, adjust_hsv);
  }

private:
  template <class S_>
  static Color3f getColor(double value, const Histogram<S_> &histogram,
                          Scheme scheme=RAINBOW, bool adjust_range=false,
                          bool adjust_hsv=false) {
    double s;
    Color3f hsv, rgb;
    
    if(adjust_range)
      s = histogram.cumulativeFrequency(value);
    else
      s = (value - histogram.min()) / histogram.rangeLength();
    
    switch(scheme) {
    case BWR:
      rgb = bwrMap(s);
      break;
    case RAINBOW:
      rgb = rainbowMap(s);
      break;
    case HEAT:
      rgb = heatMap(s);
      break;
    case GRAY:
      rgb = grayMap(s);
      break;
    }
    
    if(adjust_hsv) {
      hsv = rgbToHsv(rgb);
      s = 0.60 + (0.40 * (value - histogram.min()) / histogram.rangeLength());
      hsv[2] = s;
      rgb = hsvToRgb(hsv);
    }

    return rgb;
  }
};

#endif /* __COLOR_SCHEME_HPP__ */
