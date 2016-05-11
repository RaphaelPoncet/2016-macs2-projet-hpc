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

#include <cassert>
#include <fstream>
#include <iostream>

#include "common.hpp"
#include "location_output.hpp"
#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"

static void ValidateIndex(int i) {

  assert(0 <= i);

}

LocationOutput::LocationOutput(int index_slow):
  m_index_slow(index_slow) {

  ValidateIndex(m_index_slow);

}

void LocationOutput::WriteHeader(const RectilinearGrid3D& propagation_grid, 
                                 int nb_timesteps, std::ostream* os_ptr) {

  if (!*os_ptr) {

    LOG_ERROR << "Invalid output stream";
    std::abort();

  }

  // The variable to dump, and its type, are currently hardcoded.
  const std::string header_header = "#";
  const std::string variable_name = "pressure_0";
  const std::string variable_datatype = (sizeof(RealT) == 8 ? "Float64" : "Float32");
  const int nb_components = 1;
  const int n_fast = propagation_grid.n_fast();
  const int n_medium = propagation_grid.n_medium();
  const int n_slow = nb_timesteps;

  *os_ptr << header_header << " " 
          << variable_name << " "
          << variable_datatype << " "
          << nb_components << " "
          << n_fast << " "
          << n_medium << " "
          << n_slow << "\n";

}

void LocationOutput::Write(const RectilinearGrid3D& propagation_grid, 
                           const MultiDimensionalStorage4D& variable_storage,
                           std::ostream* os_ptr) {

  if (!*os_ptr) {

    LOG_ERROR << "Invalid output stream";
    std::abort();

  }

  assert(0);

}
