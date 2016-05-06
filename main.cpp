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
#include <sstream>

#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"
#include "variable_definitions.hpp"

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

  // Variables information (hardcoded for now, in variable_definitions.hpp).
  const int nb_variables = variable::NB_VARIABLES;
  
  std::stringstream var_name_msg;

  var_name_msg << "List of variables: [";
  for (int ivar = 0; ivar < nb_variables; ++ivar)
    var_name_msg << " " << variable::VARIABLE_NAMES[ivar];
  var_name_msg << " ]\n";

  LOG_DEBUG << var_name_msg.str();

  const int padding = 32;

  // Variables size and extent in all 3 dimensions.
  const int nx = 100;
  const RealT x_min = 0.0;
  const RealT x_max = 10000.0;

  const int ny = 1;
  const RealT y_min = 0.0;
  const RealT y_max = 0.0;

  const int nz = 100;
  const RealT z_min = 0.0;
  const RealT z_max = 10000.0;

  assert(0 < nx);
  const RealT dx = (x_max - x_min) / (RealT)nx;

  assert(0 < ny);
  const RealT dy = (y_max - y_min) / (RealT)ny;

  assert(0 < nz);
  const RealT dz = (z_max - z_min) / (RealT)nz;
  
  // Initialize grid from variable informations. 

  // We can choose our variables on grid nodes or grid
  // cells. Depending on that, grid extent and size will change a little bit.
  //
  // In our solver, all variables have the same support, so it is
  // merely a choice of convenience.

  const VariableSupport variable_support = variable::VARIABLE_SUPPORT;

  // The number of nodes is the number of cells + 1, except if the
  // number of nodes is 1, when the number of cells is also 1.
  const int nx_grid = 
    (variable_support == CELL ? (nx == 1 ? nx : nx + 1) : nx);
  const int ny_grid = 
    (variable_support == CELL ? (ny == 1 ? ny : ny + 1) : ny);
  const int nz_grid = 
    (variable_support == CELL ? (nz == 1 ? nz : nz + 1) : nz);

  const RealT x_min_grid = (variable_support == CELL ? x_min - 0.5 * dx: x_min);
  const RealT x_max_grid = (variable_support == CELL ? x_max + 0.5 * dx: x_max);
  const RealT y_min_grid = (variable_support == CELL ? y_min - 0.5 * dy: y_min);
  const RealT y_max_grid = (variable_support == CELL ? y_max + 0.5 * dy: y_max);
  const RealT z_min_grid = (variable_support == CELL ? z_min - 0.5 * dz: z_min);
  const RealT z_max_grid = (variable_support == CELL ? z_max + 0.5 * dz: z_max);

  RectilinearGrid3D propagation_grid = RectilinearGrid3D(x_min_grid, x_max_grid, nx_grid, 
							 y_min_grid, y_max_grid, ny_grid, 
							 z_min_grid, z_max_grid, nz_grid);

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
