/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __HISTOGRAM_HPP__
#define __HISTOGRAM_HPP__

#include <assert.h>

#include <iostream>
#include <limits>
#include <vector>

#include "utils/general.hpp"

template <class T_>
struct Interval {
  T_ lo, hi;

  Interval() {
    lo = hi = (T_)0;
  }

  T_ length() const {
    return hi - lo;
  }
};

template <class T_>
class Histogram {
public:
  Histogram(int num_bins) : frequencies(num_bins), cumulative(num_bins) {
    numBins = num_bins;
  }

  T_ min() const {
    return range.lo;
  }

  T_ max() const {
    return range.hi;
  }

  T_ rangeLength() const {
    return range.length();
  }

  // Build the histogram
  void build(const std::vector<T_> &data) {
    // Initialize frequencies
    for(int i = 0; i < numBins; i++) {
      frequencies[i] = 0.0f;
    }

    int num_items = data.size();
    range.lo = std::numeric_limits<T_>::max();
    range.hi = -std::numeric_limits<T_>::min();

    // Find minimum and maximum, and validate_data
    for(unsigned int i = 0; i < data.size(); i++) {
      assert(is_valid(data[i]));
      if(data[i] < range.lo) range.lo = data[i];
      if(data[i] > range.hi) range.hi = data[i];
    }

    // Expand the range a little bit
    T_ expansion = (1.0 + range.length() / numBins) * 0.01;
    range.lo -= expansion;
    range.hi += expansion;

    T_ range_length = range.length();
    T_ interval_length = range_length / numBins;

    // Count the frequencies
    for(unsigned int i = 0; i < data.size(); i++) {
      int bin = (data[i] - range.lo) / interval_length;
      frequencies[bin]++;
    }

    // Normalize and accumulate the frequencies
    frequencies[0] /= num_items;
    cumulative[0] = frequencies[0];

    for(int i = 1; i < numBins; i++) {
      frequencies[i] /= num_items;
      cumulative[i] = cumulative[i - 1] + frequencies[i];
    }
  }

  // Returns a normalized frequency for the bins
  // in which value falls
  T_ frequency(T_ value) const {
    T_ interval_length = range.length() / numBins;
    int bin = (value - range.lo) / interval_length;
    return (T_)frequencies[bin];
  }

  T_ cumulativeFrequency(T_ value) const {
    T_ interval_length = range.length() / numBins;
    int bin = (value - range.lo) / interval_length;
    T_ y_hi = cumulative[bin];
    T_ y_lo = (bin == 0 ? 0.0f : cumulative[bin - 1]);
    T_ x_hi = range.lo + interval_length * (bin + 1);
    T_ x_lo = x_hi - interval_length;

    return y_lo + (value - x_lo) * (y_hi - y_lo) / (x_hi - x_lo);
  }

private:
  Interval<T_> range;
  std::vector<T_> frequencies;
  std::vector<T_> cumulative;
  int numBins;
};

#endif
