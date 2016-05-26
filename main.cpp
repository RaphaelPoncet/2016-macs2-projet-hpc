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
#include "timeloop_manager.hpp"
#include "variable_definitions.hpp"
#include "wave_propagation.hpp"
#include "wave_solver_options.hpp"

struct timespec timer_start;
struct timespec timer_stop;

std::vector<RealT> timings_numerics;

static int Usage(const std::string& exe_name) {

  std::cout << "Usage: " << exe_name << " [parameter file]" << std::endl;
  return 0;

}

int main(int argc, char** argv) {

  if (argc != 2)
    return Usage(std::string(argv[0]));

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

  char const* const parameter_filename_c = argv[1];
  const std::string parameter_filename = std::string(parameter_filename_c);

  RectilinearGrid3D propagation_grid;
  MultiDimensionalStorage4D variable_storage;
  MathematicalParser math_parser;
  std::vector<OutputEvent> output_events;
  TimeloopManager timeloop_manager;
  WaveSolverOptions wave_solver_options = WaveSolverOptions();;

  LOG_INFO << "Reading parameter file \"" << parameter_filename << "\"...";

  std::ifstream parameter_file(parameter_filename.c_str(), std::ifstream::in);
  const int nx_padding = 0;

  ParseParameterFile(nx_padding, &parameter_file, &propagation_grid, 
                     &variable_storage, &math_parser, &output_events,
                     &timeloop_manager, &wave_solver_options);

  parameter_file.close();

  math_parser.PrintConstants();

  const int nb_iter = timeloop_manager.nb_iter();

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

  RealT* sponge_fast = VectorRawData<RealT>(wave_solver_options.sponge_fast());
  RealT* sponge_medium = VectorRawData<RealT>(wave_solver_options.sponge_medium());
  RealT* sponge_slow = VectorRawData<RealT>(wave_solver_options.sponge_slow());

  std::ofstream sponge_fast_out("sponge_fast.txt", std::ofstream::out);
  std::ofstream sponge_medium_out("sponge_medium.txt", std::ofstream::out);
  std::ofstream sponge_slow_out("sponge_slow.txt", std::ofstream::out);

  wave_solver_options.DumpAllSponges(&sponge_fast_out,
                                     &sponge_medium_out,
                                     &sponge_slow_out);

  sponge_fast_out.close();
  sponge_medium_out.close();
  sponge_slow_out.close();

  const RealT dt = timeloop_manager.dt();
  RealT t = 0.0;

  // Main time loop. All numerics (except time step computation) is in
  // here.
  for (int iter = 0; iter < timeloop_manager.nb_iter(); ++iter) {

    t += dt;
    math_parser.AddConstant("t", t);

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

    timer_stop = Timer::Now();
    
    timings_numerics.push_back(Timer::DiffTime(timer_start, timer_stop));

    for (auto it = output_events.begin(); it != output_events.end(); ++it) {

      const bool event_happens = 
        ((iter % it->rhythm() == 0) && (it->rhythm() >= 0)) ||
        ((iter == timeloop_manager.nb_iter() - 1) && (it->rhythm() == - 1));

      if (event_happens) {

        it->Execute(iter, propagation_grid, &math_parser, &variable_storage);

      }
    }

    // Exchange pressure arrays.
    std::swap(pressure_0, pressure_1);

  }

  // Variables cleanup.
  variable_storage.DeAllocate();

  for (auto it = output_events.begin(); it != output_events.end(); ++it)
    it->Destroy();

  const RealT minimum_time = 
    *std::min_element(timings_numerics.begin(), timings_numerics.end());
  assert(0.0 < minimum_time);

  const size_t nb_grid_points = 
    propagation_grid.n_fast() * propagation_grid.n_medium() * propagation_grid.n_slow();

  LOG_INFO << "****************************************";
  std::cout << "\n";
  Timer::PrintTimings(timings_numerics, "            Wave propagation", &std::cout);
  std::cout << "            Wave propagation : " 
            << std::scientific << (RealT)nb_grid_points / (minimum_time * 1.0e-3)
            << " grid updates per second";
  std::cout << "\n";
  std::cout << "\n";
  LOG_INFO << "****************************************";

  return 0;
}
