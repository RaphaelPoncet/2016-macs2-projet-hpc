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
#include "mathematical_parser.hpp"
#include "multidimensional_storage.hpp"
#include "parameter_parser.hpp"
#include "rectilinear_grid.hpp"
#include "timer.hpp"
#include "variable_definitions.hpp"
#include "wave_propagation.hpp"

static struct timespec timer_start;
static struct timespec timer_stop;
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

  RectilinearGrid3D propagation_grid = RectilinearGrid3D();
  MultiDimensionalStorage4D variable_storage;
  MathematicalParser math_parser;
  std::vector<OutputEvent> output_events;

  const std::string parameter_filename = "./homogeneous.json";
  LOG_INFO << "Reading parameter file \"" << parameter_filename << "\"...";

  std::ifstream parameter_file(parameter_filename.c_str(), std::ifstream::in);
  const int nx_padding = 17;
  ParseParameterFile(nx_padding, &parameter_file, &propagation_grid, 
                     &variable_storage, &math_parser, &output_events);
  parameter_file.close();

  math_parser.PrintConstants();

  const int nb_iter = 1000;

  // Variables information (hardcoded for now, in variable_definitions.hpp).
  const int nb_variables = variable::NB_VARIABLES;
  
  std::stringstream var_name_msg;

  var_name_msg << "List of variables: [";
  for (int ivar = 0; ivar < nb_variables; ++ivar)
    var_name_msg << " " << variable::VARIABLE_NAMES[ivar];
  var_name_msg << " ]\n";

  LOG_DEBUG << var_name_msg.str();

  const RealT dx = propagation_grid.dx_fast();
  const RealT dy = propagation_grid.dx_medium();
  const RealT dz = propagation_grid.dx_slow();

  RealT* pressure_0 = variable_storage.RawDataSlowDimension(variable::PRESSURE_0);
  RealT* pressure_1 = variable_storage.RawDataSlowDimension(variable::PRESSURE_1);
  RealT* laplace_p = variable_storage.RawDataSlowDimension(variable::LAPLACE_PRESSURE);
  RealT* velocity = variable_storage.RawDataSlowDimension(variable::VELOCITY);

  // Init output events.
  for (auto it = output_events.begin(); it != output_events.end(); ++it)
    it->Create(nb_iter, variable_storage, propagation_grid);

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

  const RealT dt = 0.0005;

  // Main time loop. All numerics (except time step computation) is in
  // here.
  for (int iter = 0; iter < nb_iter; ++iter) {

    const int stencil_radius = 1;

    const int n_fast = variable_storage.n_fast();
    const int n_fast_padding = variable_storage.n_fast_padding();
    const int n_medium = variable_storage.n2();
    const int n_slow = variable_storage.n3();

    timer_start = Timer::Now();

    // Apply sponge.
    ApplySpongeLayer_0(n_fast, n_fast_padding, n_medium, n_slow,
                       sponge_fast, sponge_medium, sponge_slow, 
                       pressure_1);

    ComputeLaplacian_0(n_fast, n_fast_padding, n_medium, n_slow,
                       stencil_radius, dx, dy, dz, 
                       pressure_0, laplace_p);
     
    // Advance pressure.
    AdvanceWavePressure_0(n_fast, n_fast_padding, n_medium, n_slow, 
                          stencil_radius, dt,
                          velocity, laplace_p, pressure_0, pressure_1);

    // Apply sponge.
    ApplySpongeLayer_0(n_fast, n_fast_padding, n_medium, n_slow,
                       sponge_fast, sponge_medium, sponge_slow, 
                       pressure_1);

    // Exchange pressure arrays.
    std::swap(pressure_0, pressure_1);

    timer_stop = Timer::Now();

    timings_update.push_back(Timer::DiffTime(timer_start, timer_stop));

    for (auto it = output_events.begin(); it != output_events.end(); ++it) {

      const bool event_happens = (iter % it->rhythm() == 0);

      if (event_happens) {

        // LOG_VERBOSE << "Iteration " << iter << ", event " << it->type() << " happens";
        it->Execute(iter, variable_storage, propagation_grid);

      }
    }
  }

  free(sponge_fast);
  free(sponge_medium);
  free(sponge_slow);

  // Variables cleanup.
  variable_storage.DeAllocate();

  for (auto it = output_events.begin(); it != output_events.end(); ++it)
    it->Destroy();

  Timer::PrintTimings(timings_update, "WavePropagation", &std::cerr);
 
  return 0;
}
