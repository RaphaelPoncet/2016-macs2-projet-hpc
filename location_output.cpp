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

#include "array_binary_io.hpp"
#include "common.hpp"
#include "location_output.hpp"
#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"
#include "variable_definitions.hpp"

static void ValidateIndex(int i, const RectilinearGrid3D& grid) {

  if ((i < 0) || (i >= grid.n_slow())) {

    LOG_ERROR << "Invalid index " << i
              << " for LocationOutput, must be between 0 and "
              << grid.n_slow() << " (grid slow dimension)";

    std::abort();

  }

}

LocationOutput::LocationOutput(int index_slow, const RectilinearGrid3D& grid):
  m_index_slow(index_slow) {

  ValidateIndex(m_index_slow, grid);

}

void LocationOutput::WriteHeader(const RectilinearGrid3D& propagation_grid, 
                                 int nb_timesteps, std::ostream* os_ptr) {

  if (!*os_ptr) {

    LOG_ERROR << "Invalid output stream";
    std::abort();

  }

  // Currently hardcoded.
  const std::string variable_name = "pressure_0";
  // The I/O data type must be coherent with the one used internally in the code.
  const DataType data_type = (sizeof(RealT) == 4 ? FLOAT32 : FLOAT64);
  const int nb_components = 1;

  // N.B. Taking for the slow dimensions of the variable the same as
  // the grid is just a convention (because the thirs variable for
  // LocationOutput is time).
  WriteBinaryVariableHeader(variable_name, data_type, 
                            GRIDDED_BINARY,
                            nb_components, 
                            propagation_grid.n_fast(), 
                            propagation_grid.n_medium(), 
                            nb_timesteps, 
                            propagation_grid.x_fast_min(), propagation_grid.x_fast_max(),
                            propagation_grid.x_medium_min(), propagation_grid.x_medium_max(),
                            propagation_grid.x_slow_min(), propagation_grid.x_slow_max(),
                            os_ptr);

}

void LocationOutput::Write(const RectilinearGrid3D& propagation_grid, 
                           const MultiDimensionalStorage4D& variable_storage,
                           std::ostream* os_ptr) {

  UNUSED(variable::VARIABLE_NAMES);
  UNUSED(propagation_grid);

  if (!*os_ptr) {

    LOG_ERROR << "Invalid output stream";
    std::abort();

  }

  const DataType data_type = (sizeof(RealT) == 4 ? FLOAT32 : FLOAT64);
  
  assert(variable_storage.dimensions().size() == 4);

  const int n_fast = variable_storage.dimensions().at(0);
  const int n_fast_padding = variable_storage.padding_fast();
  const int n2 = variable_storage.dimensions().at(1);
  const int index_slow = m_index_slow;
  
  const int variable_id = variable::PRESSURE_0;
  const RealT* data = variable_storage.RawDataSlowDimension(variable_id);
  assert(data != NULL);
  
  WriteBinaryVariableSliceSlowDimension(data_type,
                                        GRIDDED_BINARY,
                                        n_fast, n_fast_padding, n2, index_slow,
                                        data, os_ptr);

}
