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

#ifndef LOCATION_OUTPUT_HPP
#define LOCATION_OUTPUT_HPP

#include <fstream>
#include <iosfwd>

#include "multidimensional_storage.hpp"

class RectilinearGrid3D;

class LocationOutput {
public:
LocationOutput(int index_slow, const RectilinearGrid3D& grid);
  void WriteHeader(const RectilinearGrid3D& propagation_grid, 
                   int nb_timesteps, std::ostream* os_ptr);
  void Write(const RectilinearGrid3D& propagation_grid, 
             const MultiDimensionalStorage4D& variable_storage,
             std::ostream* os_ptr);
private:
  // We only support dumping the full hyperplane islow = index_slow.
  // In particular, it implies that the points where we can record the
  // variables are grid points.
  int m_index_slow;
};

#endif // LOCATION_OUTPUT_HPP
