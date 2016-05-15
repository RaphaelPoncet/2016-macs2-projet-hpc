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
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#include "common.hpp"
#include "location_output.hpp"
#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"
#include "variable_definitions.hpp"

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

  const int n_fast = variable_storage.n_fast();
  const int n_fast_padding = variable_storage.n_fast_padding();
  const int n_fast_padded = n_fast + n_fast_padding;
  const int n2 = variable_storage.n2();
  const int n3 = variable_storage.n3();
  const int index_slow = m_index_slow;
  
  const int variable_id = variable::PRESSURE_0;
  // LOG_INFO << "variable id = " << variable_id;
  const RealT* data = variable_storage.RawDataSlowDimension(variable_id);
  assert(data != NULL);

  const size_t index_base = index_slow * n2 * n_fast_padded;
  size_t index = index_base;

  RealT min = std::numeric_limits<RealT>::max();
  RealT max = std::numeric_limits<RealT>::min();

  for (int i2 = 0; i2 < n2; ++i2) {

    os_ptr->write(reinterpret_cast<const char*>(&(data[index])), n_fast * sizeof(RealT));
    
    index += n_fast_padded;
  }

}
