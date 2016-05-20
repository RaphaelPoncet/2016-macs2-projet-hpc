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
#include <vector>

#include "array_binary_io.hpp"
#include "grid_output.hpp"
#include "location_output.hpp"
#include "multidimensional_storage.hpp"
#include "rectilinear_grid.hpp"
#include "timer.hpp"
#include "variable_definitions.hpp"
#include "wave_propagation.hpp"

struct timespec timer_start;
struct timespec timer_stop;
std::vector<double> timings_update;

int main(int argc, char** argv) {

  UNUSED(argc);
  UNUSED(argv);

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

  const int nx_padding = 17;

  // Variables size and extent in all 3 dimensions.
  // const int nx = 384;
  // const RealT x_min = 0.0;
  // const RealT x_max = 9192.0;

  // const int ny = 1;
  // const RealT y_min = 0.0;
  // const RealT y_max = 0.0;

  // const int nz = 122;
  // const RealT z_min = 0.0;
  // const RealT z_max = 2904.0;

  const int nx = 1700;
  const RealT x_min = 0.0;
  const RealT x_max = 13600;

  const int ny = 1;
  const RealT y_min = 0.0;
  const RealT y_max = 0.0;

  const int nz = 350;
  const RealT z_min = 0.0;
  const RealT z_max = 2904.0;

  // const int nx = 3000;
  // const RealT x_min = 0.0;
  // const RealT x_max = 13600;

  // const int ny = 1;
  // const RealT y_min = 0.0;
  // const RealT y_max = 0.0;

  // const int nz = 3000;
  // const RealT z_min = 0.0;
  // const RealT z_max = 2904.0;


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
    MultiDimensionalStorage4D(nx, ny, nz, nb_variables, nx_padding);

  variable_storage.Allocate();

  RealT* pressure_0 = variable_storage.RawDataSlowDimension(variable::PRESSURE_0);
  RealT* pressure_1 = variable_storage.RawDataSlowDimension(variable::PRESSURE_1);
  RealT* velocity = variable_storage.RawDataSlowDimension(variable::VELOCITY);

  const RealT xmid = 0.5 * (x_min + x_max);
  const RealT ymid = 0.5 * (y_min + y_max);
  const RealT zmid = 500.0;//0.5 * (z_min + z_max);

  for (int iz = 0; iz < nz; ++iz) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int ix = 0; ix < nx; ++ix) {

        const int nx_pad = nx + nx_padding;

        const size_t index = ny * nx_pad * iz + nx_pad * iy + ix;

        const RealT x = (RealT)ix * dx;
        const RealT y = (RealT)iy * dy;
        const RealT z = (RealT)iz * dz;

        const RealT distance = 
          (x - xmid) * (x - xmid) + (y - ymid) * (y - ymid) + (z - zmid) * (z - zmid);

        pressure_0[index] = expf( - 0.003 * distance);
        pressure_1[index] = pressure_0[index];
        //velocity[index] = 1.0;
	
      }
    }
  }

  // Read velocity from a file.
  const std::string input_filename = "data/marm2_downsampled_8.xyz";
  // const std::string input_filename = "data/marmousi_smooth.xyz";

  LOG_INFO << "Reading input file \"" << input_filename << "\"...";

  std::ifstream velocity_file(input_filename.c_str(), std::ifstream::in | std::ifstream::binary);
  ReadBinaryVariable(velocity_file, nx, nx_padding, ny, nz, 1, velocity);

  LOG_INFO << "Reading input file done.\n";

  variable_storage.Validate();

  // Init sponge layers.
  const int sponge_width = 40;
  RealT* sponge_fast = (RealT*) malloc(propagation_grid.n_fast() * sizeof(RealT));
  RealT* sponge_medium = (RealT*) malloc(propagation_grid.n_medium() * sizeof(RealT));
  RealT* sponge_slow = (RealT*) malloc(propagation_grid.n_slow() * sizeof(RealT));

  InitSpongeArray(sponge_width, propagation_grid.n_fast(), sponge_fast);
  InitSpongeArray(sponge_width, propagation_grid.n_medium(), sponge_medium);
  InitSpongeArray(sponge_width, propagation_grid.n_slow(), sponge_slow);

  std::ofstream sponge_fast_out("sponge_fast.txt", std::ofstream::out);
  DumpSpongeArray(propagation_grid.n_fast(), sponge_fast, &sponge_fast_out);
  sponge_fast_out.close();

  std::ofstream sponge_medium_out("sponge_medium.txt", std::ofstream::out);
  DumpSpongeArray(propagation_grid.n_medium(), sponge_medium, &sponge_medium_out);
  sponge_medium_out.close();

  std::ofstream sponge_slow_out("sponge_slow.txt", std::ofstream::out);
  DumpSpongeArray(propagation_grid.n_slow(), sponge_slow, &sponge_slow_out);
  sponge_slow_out.close();

  const int nb_iter = 10000;
  
  const int index_slow = 100;

  LocationOutput location_output = LocationOutput(index_slow);
  const std::string receiver_filename = "output/receivers.txt";
  // For now, write receiver file every timestep.
  const int receiver_output_rhythm = 10;
  std::ofstream receiver_file;
  receiver_file.open(receiver_filename.c_str(), std::ios::out | std::ios::binary);
  location_output.WriteHeader(propagation_grid, nb_iter / receiver_output_rhythm, &receiver_file);

  const int validation_rhythm = 2000;

  const int output_rhythm = 2000;
  const std::string base_name = "output/output";
  const std::string vtk_extension = ".vtr";
  const std::string binary_extension = ".xyz";

  const RealT dt = 0.0005;

  for (int iter = 0; iter < nb_iter; ++iter) {

    const RealT one = 1.0;
    const RealT two = 2.0;
    const RealT hx = one / 24.0;
    const RealT hy = one / 24.0;
    UNUSED(hy);
    const RealT hz = one / 24.0;

    const int radius = 1;

    const int n_fast = propagation_grid.n_fast();    
    const int n_fast_pad = n_fast + variable_storage.n_fast_padding();
    const int n_fast_min = radius;
    const int n_fast_max = n_fast - radius;

    const int n_medium = propagation_grid.n_medium();
    const int n_medium_min = (n_medium == 1 ? 0 : radius);
    const int n_medium_max = (n_medium == 1 ? n_medium : n_medium - radius);

    const int n_slow = propagation_grid.n_slow();
    const int n_slow_min = radius;
    const int n_slow_max = n_slow - radius;

    // // Sponge layer.
    // for (int islow = 0; islow < n_slow; ++islow) {
    //   for (int imedium = 0; imedium < n_medium; ++imedium) {
    //     for (int ifast = 0; ifast < n_fast; ++ifast) {
          
    //       const size_t index = 
    //         n_medium * n_fast_pad * islow + n_fast_pad * imedium + ifast;

    //       const RealT scaling = sponge_slow[islow] * sponge_medium[imedium] * sponge_fast[ifast];
    //       // const RealT scaling = 1.0;

    //       pressure_1[index] *= scaling;

    //     }
    //   }
    // }

    timer_start = Timer::Now();
     
    // Advance pressure.
    for (int islow = n_slow_min; islow < n_slow_max; ++islow) {
      for (int imedium = n_medium_min; imedium < n_medium_max; ++imedium) {
        for (int ifast = n_fast_min; ifast < n_fast_max; ++ifast) {
            
          const size_t index = 
            n_medium * n_fast_pad * islow + n_fast_pad * imedium + ifast;

          const RealT s = velocity[index] * velocity[index] * dt * dt;
          // const RealT s = 1000.0 * 1000.0 * dt * dt;

          pressure_1[index] = two * pressure_0[index] - pressure_1[index];
          pressure_1[index] += s * hx * hx * (pressure_0[index + 1] - two * pressure_0[index] + pressure_0[index - 1]);
          pressure_1[index] += s * hz * hz * (pressure_0[index + n_fast_pad] - two * pressure_0[index] + pressure_0[index - n_fast_pad]);
        
        }
      }
    }

    timer_stop = Timer::Now();

    timings_update.push_back(Timer::DiffTime(timer_start, timer_stop));

    // // Sponge layer.
    // for (int islow = 0; islow < n_slow; ++islow) {
    //   for (int imedium = 0; imedium < n_medium; ++imedium) {
    //     for (int ifast = 0; ifast < n_fast; ++ifast) {
          
    //       const size_t index = 
    //         n_medium * n_fast_pad * islow + n_fast_pad * imedium + ifast;

    //       const RealT scaling = sponge_slow[islow] * sponge_medium[imedium] * sponge_fast[ifast];
    //       // const RealT scaling = 1.0;

    //       pressure_1[index] *= scaling;

    //     }
    //   }
    // }
    
    std::swap(pressure_0, pressure_1);

    if (iter % receiver_output_rhythm == 0) {
      
      location_output.Write(propagation_grid, variable_storage, &receiver_file);

    }

    if (iter % output_rhythm == 0) {

      std::stringstream vtk_sstr_output_filename;
      vtk_sstr_output_filename << base_name;
      vtk_sstr_output_filename << std::setfill('0') << std::setw(5) << iter;
      vtk_sstr_output_filename << vtk_extension;

      const std::string vtk_output_filename = vtk_sstr_output_filename.str();

      LOG_INFO << "Writing output file \"" << vtk_output_filename << "\"...";
    
      std::ofstream vtk_output_file(vtk_output_filename.c_str(), std::ofstream::binary);
      
      OutputGridAndData(propagation_grid, variable_storage, &vtk_output_file);
      
      vtk_output_file.close();
      
      LOG_INFO << "Writing output file done.\n";

      std::stringstream binary_sstr_output_filename;
      binary_sstr_output_filename << base_name;
      binary_sstr_output_filename << std::setfill('0') << std::setw(5) << iter;
      binary_sstr_output_filename << binary_extension;

      const std::string binary_output_filename = binary_sstr_output_filename.str();

      LOG_INFO << "Writing output file \"" << binary_output_filename << "\"...";
    
      std::ofstream binary_output_file(binary_output_filename.c_str(), std::ofstream::binary);
      
      WriteBinaryVariable(variable::VARIABLE_NAMES[variable::PRESSURE_0],
                          variable_storage.n_fast(), variable_storage.n_fast_padding(),
                          variable_storage.n2(), variable_storage.n3(), 1,
                          variable_storage.RawDataSlowDimension(variable::PRESSURE_0),
                          &binary_output_file);
      
      binary_output_file.close();
      
      LOG_INFO << "Writing output file done.\n";

    }



    if (iter % validation_rhythm == 0) {

      LOG_INFO << "Iteration " << iter;

      // For debug.
      variable_storage.Validate();

    }
  }

  receiver_file.close();

  free(sponge_fast);
  free(sponge_medium);
  free(sponge_slow);

  // Variables cleanup.
  variable_storage.DeAllocate();

  // Convert receiver file to VTK.

  // Time is in milliseconds.
  const RectilinearGrid3D receiver_grid = 
    RectilinearGrid3D(x_min_grid, x_max_grid, nx, 
                      y_min_grid, y_max_grid, ny,
                      z_min_grid, z_max_grid,
                      nb_iter / receiver_output_rhythm);

  MultiDimensionalStorage4D variable_storage_receivers = 
    MultiDimensionalStorage4D(nx, ny, nb_iter / receiver_output_rhythm, nb_variables, nx_padding);

  variable_storage_receivers.Allocate();
  
  RealT* data = variable_storage_receivers.RawDataSlowDimension(variable::PRESSURE_0);

  std::ifstream receiver_file_binary;
  LOG_INFO << "Reading receiver file "
           << "\'" << receiver_filename << "\'"
           << " for VTK conversion...";

  receiver_file_binary.open(receiver_filename.c_str(), std::ifstream::in | std::ifstream::binary);
  ReadBinaryVariable(receiver_file_binary, nx, nx_padding, ny, nb_iter / receiver_output_rhythm, 1, data);
  receiver_file_binary.close();

  LOG_INFO << "Reading receiver file done.";
  
  variable_storage_receivers.Validate();

  const std::string output_filename = "output/receivers.vtr";

  LOG_INFO << "Writing VTK receiver file \"" << output_filename << "\"...";
  
  std::ofstream receivers_vtk_file(output_filename.c_str(), std::ofstream::binary);
  
  OutputGridAndData(receiver_grid, variable_storage_receivers, &receivers_vtk_file);
  
  receivers_vtk_file.close();
  
  LOG_INFO << "Writing VTK receiver file done.\n";
  
  variable_storage_receivers.DeAllocate();

  LOG_INFO << "";

  Timer::PrintTimings(timings_update, "WavePropagation", &std::cerr);
 
  return 0;
}
