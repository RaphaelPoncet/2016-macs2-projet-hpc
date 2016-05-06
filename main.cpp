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

#include "common.hpp"

#include <fstream>

#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"

int main(int argc, char** argv) {

  // Initialize logging.
  const bool log_to_file = false;

  if (log_to_file) {

    const int nb_rolling_logfiles = 2;
    const int max_logfile_size = 1000000;
    plog::init(plog::verbose, "wave_propagator.log", max_logfile_size, nb_rolling_logfiles);

  } else {

    static plog::ConsoleAppender<plog::TxtFormatter> console_appender;
    plog::init(plog::verbose, &console_appender);

  }

  LOG_INFO << "**** 2D/3D wave propagator in heterogeneous media ****\n";

  // Variables Initialization (hardcoded for now).
  const int nx = 128;
  const int ny = 128;
  const int nz = 128;
  const int nb_variables = 4;
  const int padding = 32;

  const RealT x_min = 0.0;
  const RealT x_max = 10000.0;
  const RealT y_min = 0.0;
  const RealT y_max = 10000.0;
  const RealT z_min = 0.0;
  const RealT z_max = 10000.0;
  
  RectilinearGrid3D propagation_grid =
    RectilinearGrid3D(x_min, x_max, nx, y_min, y_max, ny, z_min, z_max, nz);

  std::ofstream grid_file;
  const std::string grid_filename = "grid.vtr";

  grid_file.open(grid_filename.c_str());

  propagation_grid.WriteHeaderVTKXml(&grid_file);
  propagation_grid.WriteVTKXmlAscii(&grid_file);
  propagation_grid.WriteFooterVTKXml(&grid_file);

  grid_file.close();
  
  MultiDimensionalStorage4D variable_storage = 
    MultiDimensionalStorage4D(nx, ny, nz, nb_variables, padding);

  variable_storage.Allocate();

  for (int ivar = 0; ivar < nb_variables; ++ivar) {
    for (int iz = 0; iz < nz; ++iz) {
      for (int iy = 0; iy < ny; ++iy) {
	for (int ix = 0; ix < nx; ++ix) {

	  const int nx_pad = nx + padding;

	  const size_t index = 
	    nz * ny * nx_pad * ivar + ny * nx_pad * iz + nx_pad * iy + ix;
	  variable_storage.m_data[index] = 1.0;

	}
      }
    }
  }

  variable_storage.Validate();

  // Variables cleanup.
  variable_storage.DeAllocate();

  return 0;
}
