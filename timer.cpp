// Copyright 2016 Raphael Poncet.

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//      http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language goveroning permissions and
// limitations under the License.

#include <algorithm>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "timer.hpp"

namespace Timer {

  struct timespec Now() {

    struct timespec temp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &temp);

    return temp;

  }

  double DiffTime(const struct timespec& start, const struct timespec& end) {
    
    struct timespec temp;

    if ((end.tv_nsec - start.tv_nsec) < 0) {

      temp.tv_sec = end.tv_sec-start.tv_sec - 1;
      temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;

    } else {

      temp.tv_sec = end.tv_sec-start.tv_sec;
      temp.tv_nsec = end.tv_nsec-start.tv_nsec;

    }

    const double time_in_ms = 1.0e3 * temp.tv_sec + 1.0e-6 * temp.tv_nsec;

    return time_in_ms;

  }

  void PrintTimings(const std::vector<RealT>& timings, 
                    const std::string& kernel_name, 
                    std::ostream* os_ptr) {

    if (!timings.empty()) {

      const double mean = std::accumulate(timings.begin(), timings.end(), 0.0);
  
      *os_ptr << std::setprecision(10) << kernel_name << " : "
              << "min=" 
              << *std::min_element(timings.begin(), timings.end()) << "ms, "
              << "mean=" 
              << mean / timings.size() << "ms, "
              << "max=" 
              << *std::max_element(timings.begin(), timings.end()) << "ms\n";

    }
  
  }

}
