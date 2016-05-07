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

#include <cmath>
#include <fstream>
#include <sstream>

#include "array_vtk_io.hpp"
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

    static plog::ColorConsoleAppender<plog::TxtFormatter> console_appender;
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

  // Depending on variable support (cell or node), grid extent
  // and size will change a little.
  //
  // It is merely a matter of convenience for visualization: by
  // default, Paraview outputs cell variables with no interpolation,
  // and node variables with interpolation.

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
  
  MultiDimensionalStorage4D variable_storage = 
    MultiDimensionalStorage4D(nx, ny, nz, nb_variables, padding);

  variable_storage.Allocate();

  RealT* pressure_0 = variable_storage.RawDataSlowDimension(variable::PRESSURE_0);
  RealT* pressure_1 = variable_storage.RawDataSlowDimension(variable::PRESSURE_1);

  const RealT xmid = 0.5 * (x_min + x_max);
  const RealT ymid = 0.5 * (y_min + y_max);
  const RealT zmid = 0.5 * (z_min + z_max);

  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {

	const int nx_pad = nx + padding;

	const size_t index = ny * nx_pad * iz + nx_pad * iy + ix;

	const RealT x = (RealT)ix * dx;
	const RealT y = (RealT)iy * dy;
	const RealT z = (RealT)iz * dz;

	const RealT distance = 
	  (x - xmid) * (x - xmid) + (y - ymid) * (y - ymid) + (z - zmid) * (z - zmid);

	pressure_0[index] = expf(- 0.000005 * distance);
	pressure_1[index] = 3.14f;
	
      }
    }
  }

  variable_storage.Validate();

  std::ofstream output_file;
  const std::string output_filename = "grid.vtr";

  LOG_INFO << "Writing output file \"" << output_filename << "\"...";
 
  output_file.open(output_filename.c_str(), std::ofstream::binary);

  const DataFormat format = BINARY;

  propagation_grid.WriteHeaderVTKXml(&output_file);
  propagation_grid.WriteVTKXmlAscii(&output_file);

  if (format == ASCII) {

    output_file << "<CellData>\n";
  
    if (variable_support == CELL) {

      const int n_slow = (propagation_grid.n_slow() == 1 ? 1 : propagation_grid.n_slow() - 1);
      assert(0 < n_slow);
      const int n_medium = (propagation_grid.n_medium() == 1 ? 1 : propagation_grid.n_medium() - 1);
      assert(0 < n_medium);
      const int n_fast = (propagation_grid.n_fast() == 1 ? 1 : propagation_grid.n_fast() - 1);
      assert(0 < n_fast);

      const size_t nb_grid_elements = n_slow * n_medium * n_fast;

      const size_t nb_elements_storage_3d = 
	variable_storage.n3() * variable_storage.n2() * variable_storage.n_fast();

      assert(nb_elements_storage_3d == nb_grid_elements);

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

	const std::string variable_name = variable::VARIABLE_NAMES[ivar];
	const int nb_components = 1;

	if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

	  size_t grid_data_offset_in_bytes = 0;
	  const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);
	  RealT* data = variable_storage.RawDataSlowDimension(ivar);

	  output_file << "<DataArray ";
	  WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
				    &grid_data_offset_in_bytes, &output_file);
	  output_file << ">\n";

	  WriteVTKXmlVariable(format, nb_components, 
			      variable_storage.n_fast(), variable_storage.padding(),
			      variable_storage.n2(), variable_storage.n3(), 
			      data, &grid_data_offset_in_bytes, &output_file);

	  output_file << "</DataArray>\n";

	}
      }
    }

    output_file << "</CellData>\n";

    output_file << "<PointData>\n";

    if (variable_support == NODE) {

      const size_t nb_grid_elements = 
	propagation_grid.n_slow() * propagation_grid.n_medium() * propagation_grid.n_fast();

      const size_t nb_elements_storage_3d = 
	variable_storage.n_fast() * variable_storage.n2() * variable_storage.n3();

      assert(nb_elements_storage_3d == nb_grid_elements);

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

	const std::string variable_name = variable::VARIABLE_NAMES[ivar];
	const int nb_components = 1;

	if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

	  const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);
	  size_t grid_data_offset_in_bytes = 0;     
	  RealT* data = variable_storage.RawDataSlowDimension(ivar);
	  const int nb_components = 1;

	  output_file << "<DataArray ";
	  WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
				    &grid_data_offset_in_bytes, &output_file);
	  output_file << ">\n";

	  WriteVTKXmlVariable(format, nb_components, 
			      variable_storage.n_fast(), variable_storage.padding(),
			      variable_storage.n2(), variable_storage.n3(), 
			      data, &grid_data_offset_in_bytes, &output_file);

	  output_file << "</DataArray>\n";

	}
      }
    }

    output_file << "</PointData>\n";

  } else if (format == BINARY) {

    std::stringstream cell_data_arrays_header;
    std::stringstream point_data_arrays_header;

    cell_data_arrays_header << "<CellData>\n";
    point_data_arrays_header << "<PointData>\n";

    size_t grid_data_offset_in_bytes = 0;

    for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

      const std::string variable_name = variable::VARIABLE_NAMES[ivar];
      const int nb_components = 1;

      if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

	if (variable_support == CELL) {

	  const int n_slow = (propagation_grid.n_slow() == 1 ? 1 : propagation_grid.n_slow() - 1);
	  assert(0 < n_slow);
	  const int n_medium = (propagation_grid.n_medium() == 1 ? 1 : propagation_grid.n_medium() - 1);
	  assert(0 < n_medium);
	  const int n_fast = (propagation_grid.n_fast() == 1 ? 1 : propagation_grid.n_fast() - 1);
	  assert(0 < n_fast);

	  const size_t nb_grid_elements = n_slow * n_medium * n_fast;

	  const size_t nb_elements_storage_3d = 
	    variable_storage.n3() * variable_storage.n2() * variable_storage.n_fast();

	  assert(nb_elements_storage_3d == nb_grid_elements);

	  const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);

	  cell_data_arrays_header << "<DataArray ";
	  WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
				    &grid_data_offset_in_bytes, &cell_data_arrays_header);
	  cell_data_arrays_header << ">\n";
	  cell_data_arrays_header << "</DataArray>\n";

	} else if (variable_support == NODE) {

	  const size_t nb_grid_elements = 
	    propagation_grid.n_slow() * propagation_grid.n_medium() * propagation_grid.n_fast();

	  const size_t nb_elements_storage_3d = 
	    variable_storage.n_fast() * variable_storage.n2() * variable_storage.n3();

	  assert(nb_elements_storage_3d == nb_grid_elements);

	  const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);

	  point_data_arrays_header << "<DataArray ";
	  WriteVTKXmlVariableHeader(variable_name, format, nb_components, grid_data_size_in_bytes, 
				    &grid_data_offset_in_bytes, &point_data_arrays_header);
	  point_data_arrays_header << ">\n";
	  point_data_arrays_header << "</DataArray>\n";

	}
      }
    }

    cell_data_arrays_header << "</CellData>\n";
    point_data_arrays_header << "</PointData>\n";

    output_file << cell_data_arrays_header.str();
    output_file << point_data_arrays_header.str();
    
    output_file << "</Piece>\n"
		<< "</RectilinearGrid>\n";

    // The underscore is intentional.
    output_file << "<AppendedData encoding=\"raw\">\n_";
  
    if (variable_support == CELL) {

      const int n_slow = (propagation_grid.n_slow() == 1 ? 1 : propagation_grid.n_slow() - 1);
      assert(0 < n_slow);
      const int n_medium = (propagation_grid.n_medium() == 1 ? 1 : propagation_grid.n_medium() - 1);
      assert(0 < n_medium);
      const int n_fast = (propagation_grid.n_fast() == 1 ? 1 : propagation_grid.n_fast() - 1);
      assert(0 < n_fast);

      const size_t nb_grid_elements = n_slow * n_medium * n_fast;

      const size_t nb_elements_storage_3d = 
	variable_storage.n3() * variable_storage.n2() * variable_storage.n_fast();

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

	const int nb_components = 1;
	  
	size_t grid_data_offset_in_bytes = 0;
	const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);
	RealT* data = variable_storage.RawDataSlowDimension(ivar);

	if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

	  WriteVTKXmlVariable(format, nb_components, 
			      variable_storage.n_fast(), variable_storage.padding(),
			      variable_storage.n2(), variable_storage.n3(), 
			      data, &grid_data_offset_in_bytes, &output_file);

	}
      }

    } else if (variable_support == NODE) {

      const size_t nb_grid_elements = 
	propagation_grid.n_slow() * propagation_grid.n_medium() * propagation_grid.n_fast();

      const size_t nb_elements_storage_3d = 
	variable_storage.n_fast() * variable_storage.n2() * variable_storage.n3();

      assert(nb_elements_storage_3d == nb_grid_elements);

      for (int ivar = 0; ivar < variable::NB_VARIABLES; ++ivar) {

	const int nb_components = 1;
	  
	size_t grid_data_offset_in_bytes = 0;
	const size_t grid_data_size_in_bytes = nb_grid_elements * sizeof(RealT);
	RealT* data = variable_storage.RawDataSlowDimension(ivar);

	if (variable::VARIABLE_FLAGS[ivar] & WRITTEN) {

	  WriteVTKXmlVariable(format, nb_components, 
			      variable_storage.n_fast(), variable_storage.padding(),
			      variable_storage.n2(), variable_storage.n3(), 
			      data, &grid_data_offset_in_bytes, &output_file);

	}
      }
    }

    output_file << "\n</AppendedData>\n";

  }

  propagation_grid.WriteFooterVTKXml(&output_file);

  output_file.close();

  LOG_INFO << "Writing output file done.\n";
  
  // Variables cleanup.
  variable_storage.DeAllocate();
 
 return 0;
}
