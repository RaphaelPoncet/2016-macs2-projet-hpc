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

#ifndef PARAMETER_PARSER_HPP
#define PARAMETER_PARSE_HPP

#include <iosfwd>
#include <string>

#include "picojson.h"

#include "common.hpp"
#include "multidimensional_storage.hpp"

class RectilinearGrid3D;
class LocationOutput;
class MathematicalParser;
class TimeloopManager;
class WaveSolverOptions;

class OutputEvent {
public:
  OutputEvent();
  std::string type() { return m_type; }
  std::string type() const { return m_type; }
  std::string format() { return m_format; }
  std::string format() const { return m_format; }
  std::string ascii_or_binary() { return m_ascii_or_binary; }
  std::string ascii_or_binary() const { return m_ascii_or_binary; }
  std::string stream_name() { return m_stream_name; }
  std::string stream_name() const { return m_stream_name; }
  int rhythm() { return m_rhythm; }
  int rhythm() const { return m_rhythm; }
  std::string var_name() { return m_var_name; }
  std::string var_name() const { return m_var_name; }
  std::string formula() { return m_formula; }
  std::string formula() const { return m_formula; }
  void ParseFromJSON(const picojson::value& v);
  void Validate();
  void Execute(int iter,
               const RectilinearGrid3D& grid,
               MathematicalParser* parser_ptr,
               MultiDimensionalStorage4D* variable_storage_ptr);
  void Create(int nb_iter,
              const MultiDimensionalStorage4D& storage,
              const RectilinearGrid3D& grid);
  void Destroy();
private:
  std::string m_type;
  std::string m_format;
  std::string m_ascii_or_binary;
  std::string m_stream_name;
  int m_rhythm;
  std::ostream* m_output_stream_ptr;
  std::ofstream* m_file_ptr;
  LocationOutput* m_location_output_ptr;
  int m_index_slow;
  std::string m_var_name;
  std::string m_formula;
};

void ParseParameterFile(int n_fast_padding,
                        std::istream* input_stream_ptr,
                        RectilinearGrid3D* grid_ptr,
                        MultiDimensionalStorage4D* storage_ptr,
                        MathematicalParser* math_parser_ptr,
                        std::vector<OutputEvent>* output_events_ptr,
                        TimeloopManager* timeloop_manager_ptr,
                        WaveSolverOptions* wave_solver_options_ptr);

#endif // PARAMETER_PARSER_HPP
