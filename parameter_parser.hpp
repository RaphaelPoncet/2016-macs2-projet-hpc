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

class RectilinearGrid3D;
class MultiDimensionalStorage4D;
class LocationOutput;

class OutputEvent {
public:
  OutputEvent();
  std::string type() { return m_type; }
  std::string type() const { return m_type; }
  std::string format() { return m_format; }
  std::string format() const { return m_format; }
  std::string stream_name() { return m_stream_name; }
  std::string stream_name() const { return m_stream_name; }
  int rhythm() { return m_rhythm; }
  int rhythm() const { return m_rhythm; }
  void ParseFromJSON(const picojson::value& v);
  void Validate();
  void Execute(int iter,
               const MultiDimensionalStorage4D& variable_storage,
               const RectilinearGrid3D& grid);
  void Create(int nb_iter,
              const MultiDimensionalStorage4D& storage,
              const RectilinearGrid3D& grid);
  void Destroy();
private:
  std::string m_type;
  std::string m_format;
  std::string m_stream_name;
  int m_rhythm;
  std::ostream* m_output_stream_ptr;
  std::ofstream* m_file_ptr;
  LocationOutput* m_location_output_ptr;
  int m_index_slow;
};

void ParseParameterFile(int n_fast_padding,
                        std::istream* input_stream_ptr,
                        RectilinearGrid3D* grid_ptr,
                        MultiDimensionalStorage4D* storage_ptr,
                        std::vector<OutputEvent>* output_events_ptr);

#endif // PARAMETER_PARSER_HPP
